import subprocess
import os

class QIIME2Processor:
    def __init__(self, config):
        self.config = config
        self._set_paths()
        self._create_output_directories()

    def _set_paths(self):
        # Base paths
        self.base_output_dir = self.config['output_base_dir']
        
        # Paths for input data - should point to the main project's 'Datesets' or similar structure
        self.input_data_dir = self.config['input_data_dir']
        self.raw_fastq_dir = self.config.get('raw_fastq_dir', os.path.join(self.input_data_dir, 'raw_fastqs')) # Optional: where raw fastqs might be
        self.reference_db_dir = self.config.get('reference_db_dir', os.path.join(self.input_data_dir, 'reference_db')) # Optional: where reference dbs might be

        # Specific paths based on pool (if applicable)
        self.pool = self.config.get('pool', None) # Pool is now optional, depends on if demux is needed
        self.output_dir = os.path.join(self.base_output_dir, "qiime2_processing") # Centralized output for this script

        # General metadata and taxonomy files from the main 'Datesets' directory
        self.metadata_file = os.path.join(self.input_data_dir, "metadata.csv") # Assuming metadata.csv is at the base of input_data_dir
        self.taxonomy_tsv = os.path.join(self.input_data_dir, "taxonomy.tsv") # Assuming taxonomy.tsv is at the base of input_data_dir

        # Paths for the original wild-again structure (for reference, adjust as needed)
        self.wild_again_metadata_base_dir = self.config.get('wild_again_metadata_base_dir', 'wild-again/data/metadata')
        self.wild_again_ref_db_base_dir = self.config.get('wild_again_ref_db_base_dir', 'wild-again/data/db')
        self.wild_again_dataset_base_dir = self.config.get('wild_again_dataset_base_dir', 'wild-again/data/raw')

        # Construct specific paths based on pool if it exists
        if self.pool:
            self.pool_metadata_dir_wild_again = os.path.join(self.wild_again_metadata_base_dir, self.pool)
            self.pool_dataset_dir_wild_again = os.path.join(self.wild_again_dataset_base_dir, self.pool)
        else:
            self.pool_metadata_dir_wild_again = None
            self.pool_dataset_dir_wild_again = None


    def _create_output_directories(self):
        # Create main output directory for this processing run
        os.makedirs(self.output_dir, exist_ok=True)

        # Recreate the directory structure from the original script
        # Subdirectories that are always created
        general_sub_dirs = [
            "4-dada2", "5-taxonomy", "6-filter", "7-phylogeny", "8-export"
        ]
        for sub_dir in general_sub_dirs:
            os.makedirs(os.path.join(self.output_dir, sub_dir), exist_ok=True)

        # Create pool-specific subdirectories if a pool is specified
        if self.pool:
            pool_specific_sub_dirs = [
                "1-demux", "2-combine", "3-import"
            ]
            for sub_dir in pool_specific_sub_dirs:
                os.makedirs(os.path.join(self.output_dir, sub_dir, self.pool), exist_ok=True)

            # Specific sub-directories for merged tables within dada2
            os.makedirs(os.path.join(self.output_dir, "4-dada2", self.pool), exist_ok=True)
            
        # Create specific directories for 8-export
        os.makedirs(os.path.join(self.output_dir, "8-export", "tree"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "8-export", "rooted-tree"), exist_ok=True)


    def activate_conda_env(self):
        print(f"Activating conda environment: {self.config['conda_env']}")
        # This part is tricky to do from within a Python script directly
        # The 'conda activate' command is designed to modify the current shell's environment.
        # For subprocesses, it's generally recommended to use `conda run -n env_name command`
        # However, if QIIME commands are not found, this might be a reason.
        # For now, we add a subprocess call to 'activate' which might not persist.
        # The ultimate solution is to ensure 'qiime' is in the PATH when this script is run,
        # or prepend `conda run -n {self.config['conda_env']}` to all qiime commands.
        # For simplicity in this refactoring, we'll assume qiime is accessible.
        # self.run_command(["conda", "activate", self.config['conda_env']], description="Activating Conda Environment")
        # Skipping direct conda activate as it doesn't work well with subprocess.run directly.
        # Assuming the script is run from an already activated environment or `qiime` is in PATH.
        pass


    def run_command(self, cmd, description="Running command"):
        print(f"\n--- {description} ---")
        print(f"Command: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print("Stdout:\n", result.stdout)
            if result.stderr:
                print("Stderr:\n", result.stderr)
            return result
        except subprocess.CalledProcessError as e:
            print(f"Error during {description}:")
            print("Stdout:\n", e.stdout)
            print("Stderr:\n", e.stderr)
            raise e
        except FileNotFoundError:
            print(f"Error: Command not found. Make sure all necessary tools are installed and in your PATH.")
            raise

    def demultiplex(self):
        print("Starting demultiplexing with cutadapt...")

        if not self.pool:
            raise ValueError("POOL must be set in config for demultiplexing step.")

        # Input raw fastq files
        # The actual raw fastq file names might vary, so these should ideally be parameterized or discovered.
        # For now, matching the original script's assumption.
        forward_fastq = os.path.join(self.raw_fastq_dir, "forward.fastq.gz")
        reverse_fastq = os.path.join(self.raw_fastq_dir, "reverse.fastq.gz")

        # Input primer files - these should ideally be in self.reference_db_dir
        fw_primers = self.config.get('cutadapt_fw_primers')
        rv_primers = self.config.get('cutadapt_rv_primers')

        if not fw_primers or not os.path.exists(fw_primers):
            print(f"Error: Forward primer file not specified or found at {fw_primers}. Please update config['cutadapt_fw_primers'].")
            return
        if not rv_primers or not os.path.exists(rv_primers):
            print(f"Error: Reverse primer file not specified or found at {rv_primers}. Please update config['cutadapt_rv_primers'].")
            return
        
        if not os.path.exists(forward_fastq):
            print(f"Error: Forward fastq file not found at {forward_fastq}. Please ensure raw fastq files are in the specified directory.")
            return
        if not os.path.exists(reverse_fastq):
            print(f"Error: Reverse fastq file not found at {reverse_fastq}. Please ensure raw fastq files are in the specified directory.")
            return


        # Output directory for demultiplexing
        demux_output_dir = os.path.join(self.output_dir, "1-demux", self.pool)

        # Round 1
        cmd_round1 = [
            "cutadapt",
            "--cores", str(self.config['cutadapt_cores']),
            "--pair-adapters",
            "-g", f"file:{fw_primers}",
            "-G", f"file:{rv_primers}",
            "-o", os.path.join(demux_output_dir, "round1-{name}-1.fastq.gz"),
            "-p", os.path.join(demux_output_dir, "round1-{name}-2.fastq.gz"),
            forward_fastq, reverse_fastq,
            "--no-indels",
            "--errors", str(self.config['cutadapt_errors'])
        ]
        self.run_command(cmd_round1, "Cutadapt Round 1 Demultiplexing")

        # Round 2
        # Input for round 2 are the 'unknown' files from round 1
        # IMPORTANT: The original script had round1-unknown-2.fastq.gz then round1-unknown-1.fastq.gz
        # This seems like a swap for R1 and R2 for the unknown files.
        # I'll replicate that behavior for now.
        round1_unknown_f = os.path.join(demux_output_dir, "round1-unknown-2.fastq.gz") 
        round1_unknown_r = os.path.join(demux_output_dir, "round1-unknown-1.fastq.gz")

        # Check if unknown files exist, otherwise round 2 might fail or be unnecessary
        if not os.path.exists(round1_unknown_f) or not os.path.exists(round1_unknown_r):
            print(f"Warning: Unknown files from Cutadapt Round 1 not found in {demux_output_dir}. Skipping Cutadapt Round 2.")
            return

        cmd_round2 = [
            "cutadapt",
            "--cores", str(self.config['cutadapt_cores']),
            "--pair-adapters",
            "-g", f"file:{fw_primers}",
            "-G", f"file:{rv_primers}",
            "-o", os.path.join(demux_output_dir, "round2-{name}-1.fastq.gz"),
            "-p", os.path.join(demux_output_dir, "round2-{name}-2.fastq.gz"),
            round1_unknown_f, round1_unknown_r,
            "--no-indels",
            "--errors", str(self.config['cutadapt_errors'])
        ]
        self.run_command(cmd_round2, "Cutadapt Round 2 Demultiplexing")

    def combine_reads(self):
        print("Starting combining reads...")
        # Placeholder: This will be implemented in detail later
        pass

    def import_data(self):
        print("Starting data import into QIIME2...")
        # Placeholder: This will be implemented in detail later
        pass

    def dada2_denoising(self):
        print("Starting DADA2 denoising...")
        # Placeholder: This will be implemented in detail later
        pass

    def classify_taxonomy(self):
        print("Starting taxonomy classification...")
        # Placeholder: This will be implemented in detail later
        pass

    def filter_features(self):
        print("Starting feature filtering...")
        # Placeholder: This will be implemented in detail later
        pass

    def build_phylogeny(self):
        print("Starting phylogenetic tree construction...")
        # Placeholder: This will be implemented in detail later
        pass

    def export_data(self):
        print("Starting data export...")
        # Placeholder: This will be implemented in detail later
        pass

    def process(self):
        self.activate_conda_env()
        self.demultiplex()
        self.combine_reads()
        self.import_data()
        self.dada2_denoising()
        self.classify_taxonomy()
        self.filter_features()
        self.build_phylogeny()
        self.export_data()
        print("\nQIIME2 processing complete!")

if __name__ == "__main__":
    # Example Configuration - This needs to be adapted to the user's specific Datasets folder structure
    # The 'Datasets' folder in the user's main project directory is likely analogous to 'wild-again/data/raw' and 'wild-again/data/metadata'
    # Need to verify content of the main 'Datasets' to correctly map
    config = {
        'conda_env': 'qiime2-2023.2', # Name of the QIIME2 conda environment
        'pool': 'CB1', # Example pool, can be configured for CB1, CB2, CB3 (optional, for demultiplexing)
        'output_base_dir': 'analysis_outputs', # Relative to current working directory
        'input_data_dir': 'Datesets', # Points to the main project's 'Datesets' directory
        'raw_fastq_dir': 'Datesets/raw_fastqs', # Assuming raw fastqs are in Datesets/raw_fastqs (adjust as needed)
        'reference_db_dir': 'Datesets/reference_db', # Assuming reference DBs are in Datesets/reference_db (adjust as needed)

        # Cutadapt parameters
        'cutadapt_cores': 10,
        'cutadapt_errors': 0.05,
        'cutadapt_fw_primers': 'wild-again/data/metadata/CB3/bcpr-fw.fasta', # Placeholder, ideally in reference_db_dir
        'cutadapt_rv_primers': 'wild-again/data/metadata/CB3/bcpr-rv.fasta', # Placeholder, ideally in reference_db_dir
        'manifest_file': 'wild-again/data/metadata/CB3/manifest.csv', # Placeholder, ideally generated or in raw_fastq_dir
        'cat_script_path': 'wild-again/data/metadata/CB3/cat.sh', # Placeholder, needs to be adapted


        # DADA2 parameters
        'dada2_trim_left_f': 19,
        'dada2_trim_left_r': 0,
        'dada2_trunc_len_f': 139,
        'dada2_trunc_len_r': 212,
        'dada2_n_threads': 10,

        # Filter parameters
        'min_frequency': 1000,
        'taxa_to_exclude': ['mitochondria', 'chloroplast'],
        'taxa_to_include': ['bacteria']
    }

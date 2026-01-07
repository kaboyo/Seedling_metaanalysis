#!/bin/bash

# A script to add, commit, and push changes to a Git repository.

# Add all changes to staging
git add .

# Prompt for a commit message
echo "Enter commit message: "
read commit_message

# Commit with the provided message
git commit -m "$commit_message"

# Push to the main branch of the origin remote
# If your branch is not 'main', you can change it here.
git push origin main

echo "Changes pushed to remote repository."

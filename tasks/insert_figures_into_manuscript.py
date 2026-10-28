# usage: insert_figs.py <in.docx> <out.docx> <by_analysis dir>
# Replaces the live image above each "Figure N." / "Figure SN." legend with 00_main_figures/Figure_N.tiff
# or 00_supplementary_figures/Figure_SN.tiff, scaled to page width (max 9 in tall). Text is not touched.
import re,sys,zipfile,io
from PIL import Image
from lxml import etree
src,dst,R=sys.argv[1:4]
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'; A='http://schemas.openxmlformats.org/drawingml/2006/main'; RR='http://schemas.openxmlformats.org/officeDocument/2006/relationships'; WP='http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'
def q(t): return '{%s}%s'%(W,t)
z=zipfile.ZipFile(src); root=etree.fromstring(z.read('word/document.xml')); rels=z.read('word/_rels/document.xml.rels').decode()
def tgt(rid): m=re.search(r'<Relationship [^>]*Id="%s"[^>]*>'%rid,rels); return re.search(r'Target="([^"]+)"',m.group(0)).group(1)
def live(t): return not any(a.tag==q('del') for a in t.iterancestors())
new={}; last=None
for p in root.iter(q('p')):
    deleted=p.find(q('pPr')+'/'+q('rPr')+'/'+q('del')) is not None
    for d in p.iter('{%s}inline'%WP):
        if not deleted and d.find('.//{%s}blip'%A) is not None: last=d
    s=''.join(t.text or '' for t in p.iter(q('t')) if live(t))
    m=re.match(r'Figure (S?)(\d+)\. ',s)
    if m and last is not None:
        f='%s/%s/Figure_%s%s.tiff'%(R,'00_supplementary_figures' if m.group(1) else '00_main_figures',m.group(1),m.group(2))
        try: im=Image.open(f).convert('RGB')
        except FileNotFoundError: print('no file for',m.group(0)); last=None; continue
        w,h=im.size; buf=io.BytesIO(); im.resize((2400,round(h*2400/w)),Image.LANCZOS).save(buf,'PNG',optimize=True)
        media='word/'+tgt(last.find('.//{%s}blip'%A).get('{%s}embed'%RR)); new[media]=buf.getvalue()
        cx=6120130; cy=round(cx*h/w)
        if cy>8229600: cy=8229600; cx=round(cy*w/h)
        for e in [last.find('{%s}extent'%WP)]+list(last.iter('{%s}ext'%A)):
            if e is not None and e.get('cx'): e.set('cx',str(cx)); e.set('cy',str(cy))
        print('Figure %s%s  %.2f x %.2f in'%(m.group(1),m.group(2),cx/914400,cy/914400)); last=None
new['word/document.xml']=etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
out=zipfile.ZipFile(dst,'w',zipfile.ZIP_DEFLATED)
for i in z.infolist(): out.writestr(i,new.get(i.filename,z.read(i.filename)))
out.close(); print('replaced',len(new)-1,'images')

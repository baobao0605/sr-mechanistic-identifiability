#!/usr/bin/env python3
"""Check deliverable vector geometry, typography and source availability."""
from pathlib import Path
import json, subprocess, xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
import fitz

ROOT=Path(__file__).resolve().parents[1]
NS={'s':'http://www.w3.org/2000/svg'}
caps=json.loads((ROOT/'source_data/redrawn_figure_captions.json').read_text())
registry=json.loads((ROOT/'source_data/redrawn_figure_registry.json').read_text())

def check(r):
    p=ROOT/r['svg'];tree=ET.parse(p);svg=tree.getroot()
    texts={z.attrib['id']:''.join(z.itertext()) for z in svg.findall('s:text',NS)}
    fonts=[float(z.attrib['font-size']) for z in svg.findall('s:text',NS)]
    view=[float(x) for x in svg.attrib['viewBox'].split()]
    pxw=170/25.4*96;pxh=pxw*view[3]/view[2]
    out=subprocess.run(['inkscape',str(p),'--query-all'],check=True,capture_output=True,text=True).stdout
    bounds={}
    for line in out.splitlines():
        parts=line.split(',')
        if len(parts)==5 and parts[0] in texts:bounds[parts[0]]=list(map(float,parts[1:]))
    overflow=[];overlaps=[]
    for key,(x,y,w,h) in bounds.items():
        if min(x,y)<-.1 or x+w>pxw+.1 or y+h>pxh+.1:overflow.append({'text':texts[key],'bounds':[x,y,w,h]})
    keys=list(bounds)
    for i,k in enumerate(keys):
        x,y,w,h=bounds[k]
        for l in keys[i+1:]:
            a,b,c,d=bounds[l];iw=min(x+w,a+c)-max(x,a);ih=min(y+h,b+d)-max(y,b)
            if iw>1.0 and ih>1.0:overlaps.append({'a':texts[k],'b':texts[l],'intersection':[iw,ih]})
    pdf=fitz.open(ROOT/r['pdf']);pg=pdf[0]
    rec={**r,'svg_text_elements':len(texts),'svg_raster_images':len(svg.findall('.//s:image',NS)),'min_font_pt_measured':min(fonts)*170/680*72/25.4,'pdf_pages':len(pdf),'pdf_raster_images':len(pg.get_images(full=True)),'pdf_vector_paths':len(pg.get_drawings()),'out_of_bounds_text':overflow,'text_pair_overlaps':overlaps}
    pdf.close();return rec

with ThreadPoolExecutor(max_workers=3) as pool:result=list(pool.map(check,registry))
missing=[]
for no,cap in caps.items():
    for path in cap['sources']:
        if not (ROOT/path).exists():missing.append({'figure':no,'path':path})
report={'figures':result,'missing_caption_sources':missing,'complete_bilingual_set':len(result)==16,'all_svg_editable':all(r['svg_text_elements']>0 for r in result),'all_pdf_vector':all(r['pdf_raster_images']==0 and r['pdf_vector_paths']>0 for r in result),'all_text_within_canvas':all(not r['out_of_bounds_text'] for r in result),'all_text_pairs_clear':all(not r['text_pair_overlaps'] for r in result)}
(ROOT/'reviews').mkdir(exist_ok=True)
(ROOT/'reviews/redrawn_figure_qa.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='figures'},ensure_ascii=False))
for r in result:
    if r['out_of_bounds_text'] or r['text_pair_overlaps']:print(r['stem'],json.dumps({'overflow':r['out_of_bounds_text'],'overlap':r['text_pair_overlaps']},ensure_ascii=False))

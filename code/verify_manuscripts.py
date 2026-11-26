"""Check submission assembly, protected prose, cross-references and vector assets."""
from pathlib import Path
import re,json,hashlib,xml.etree.ElementTree as ET
import fitz

ROOT=Path(__file__).resolve().parents[1]
M=ROOT/'manuscripts'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def strip_math(s):return re.sub(r'\$\$.*?\$\$|\$(?:[^$]|\\\$)*?\$','',s,flags=re.S)

def main():
    result={'date':'2026-11-26','manuscripts':[],'figures':[],'protected_author_text':{},'all_checks_passed':False}
    caps=json.loads((ROOT/'source_data/redrawn_figure_captions.json').read_text())
    refs=(M/'references.md').read_text()
    ids=[int(x) for x in re.findall(r'(?m)^(\d+)\. ',refs)]
    assert ids==list(range(1,28)),ids
    for lang in ['ZH','EN']:
        source=(M/f'Manuscript_{lang}.md').read_text()
        core=(M/f'parts/author_core_{lang.lower()}.md').read_text().strip()
        assert source.count(core)==1,'Author core must be retained exactly once.'
        result['protected_author_text'][lang]={'verbatim':True,'sha256':hashlib.sha256(core.encode()).hexdigest()}
        assert re.findall(r'(?m)^### ([1-6])\. ',source)==list('123456')
        assert len(re.findall(r'(?m)^## (Methods|方法)$',source))==1
        assert '0009-0006-4111-2000' in source
        for i in range(1,9):assert source.count(f'<!--FIGURE_START:{i:02d}-->')==1
        theory=(ROOT/f'supplement/Theory_Appendix_{lang}.md').read_text()
        theorem_ids=re.findall(r'(?m)^### [^\n]*\b(T\d\d)\b',theory)
        assert theorem_ids==[f'T{i:02d}' for i in range(1,29)],theorem_ids
        emp=(ROOT/f'supplement/Supplementary_Empirical_{lang}.md').read_text()
        assert re.findall(r'(?m)^## E(\d+)\.',emp)==[str(i) for i in range(1,12)]
        for kind in ['Manuscript','Supplementary']:
            name=f'{kind}_{lang}';md=M/(name+'.md');pdf=M/(name+'.pdf');tex=M/(name+'.tex')
            s=md.read_text();assert len(re.findall(r'(?m)^\$\$$',s))%2==0
            # Citation ranges are checked after removing mathematical intervals.
            for c in re.findall(r'\[([1-9][0-9]*(?:[–,\-][0-9]+)*)\]',strip_math(s)):
                assert all(1<=int(x)<=27 for x in re.findall(r'\d+',c)),(name,c)
            d=fitz.open(pdf);text='\n'.join(p.get_text() for p in d)
            assert 'Yongzhang Ji' in text
            assert '(-)' not in text,(name,'Literal minus-delimiter artifact')
            log=(ROOT/'build'/(name+'.log')).read_text(errors='replace')
            assert 'Overfull' not in log and 'Missing character:' not in log
            assert 'undefined references' not in log
            assert tex.exists()
            item={'file':pdf.name,'pages':len(d),'bytes':pdf.stat().st_size,'sha256':digest(pdf),'overfull_boxes':0,'missing_glyphs':0,
                  'editable_markdown':md.name,'editable_tex':tex.name,'bookmark_count':len(d.get_toc())}
            if kind=='Manuscript':
                item['english_word_tokens']=len(re.findall(r"\b[A-Za-z][A-Za-z’'-]*\b",strip_math(s))) if lang=='EN' else None
                item['han_characters']=len(re.findall(r'[\u4e00-\u9fff]',s)) if lang=='ZH' else None
            result['manuscripts'].append(item)
    for i,item in caps.items():
        for f in item['sources']:assert (ROOT/f).is_file(),f
        for lang in ['zh','en']:
            stem=item['stem']+'_'+lang;svg=ROOT/'figures/svg'/(stem+'.svg');pdf=ROOT/'figures/pdf'/(stem+'.pdf')
            tree=ET.parse(svg);els=list(tree.getroot().iter())
            assert not any(e.tag.endswith('}image') for e in els),stem
            texts=sum(e.tag.endswith('}text') for e in els)
            doc=fitz.open(pdf);assert len(doc)==1 and not doc[0].get_images(full=True),stem
            assert (ROOT/'figures/png'/(stem+'.png')).exists()
            result['figures'].append({'figure':int(i),'language':lang,'editable_text_elements':texts,'pdf_raster_images':0,'svg_sha256':digest(svg)})
    result['numbered_references']=len(ids);result['theory_results']=28;result['empirical_sections']=11
    result['all_checks_passed']=True
    (ROOT/'source_data/submission_integrity.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['figures','protected_author_text']},ensure_ascii=False,indent=2))

if __name__=='__main__':main()

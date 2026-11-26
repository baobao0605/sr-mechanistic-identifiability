"""Assemble the bilingual manuscripts and standalone supplementary information.

Run from any directory. Per-section Markdown and figure-caption JSON are the
editable sources. The build script separately converts SVG figures to embedded
vector PDFs for typesetting.
"""
from pathlib import Path
import json,re

ROOT=Path(__file__).resolve().parents[1]
M=ROOT/'manuscripts'
P=M/'parts'

def no_yaml(s):
    return re.sub(r'\A---\n.*?\n---\n','',s,count=1,flags=re.S).strip()

def no_title(s):
    return re.sub(r'\A# [^\n]+\n','',no_yaml(s),count=1).strip()

def main():
    caps=json.loads((ROOT/'source_data/redrawn_figure_captions.json').read_text())
    refs=(M/'references.md').read_text().strip()
    for lang in ['zh','en']:
        upper=lang.upper()
        s='\n\n'.join((P/f'{part}_{lang}.md').read_text().strip()
                    for part in ['front','results_1_3','results_4_5','results_6','discussion'])
        s+='\n\n## '+('方法' if lang=='zh' else 'Methods')+'\n\n'
        s+='\n\n'.join((P/f'{part}_{lang}.md').read_text().strip()
                       for part in ['methods_math','methods_empirical','declarations'])
        s+='\n\n## '+('参考文献' if lang=='zh' else 'References')+'\n\n'+refs+'\n'
        for i in range(1,9):
            item=caps[str(i)];marker=f'<!--FIG{i:02d}-->'
            if s.count(marker)!=1:raise ValueError(f'{lang}: {marker} count {s.count(marker)}')
            label=f'图 {i}' if lang=='zh' else f'Figure {i}'
            block=(f'<!--FIGURE_START:{i:02d}-->\n\n'
                   f'![{label}](../figures/svg/{item["stem"]}_{lang}.svg)\n\n'
                   f'{item[lang]}\n\n<!--FIGURE_END:{i:02d}-->')
            s=s.replace(marker,block)
        (M/f'Manuscript_{upper}.md').write_text(s)
        preview=s.replace('../figures/svg/','../figures/png/').replace('.svg)','.png)')
        (M/f'Manuscript_{upper}_preview.md').write_text(preview)
        if lang=='zh':
            title='补充材料：认知预测地图中方向与时间权重的机制识别'
            lead='本文件提供正文各结论的一般性证明、真实重复测量与已知变换校准的完整方法和结果，以及相关研究背景。理论结果 T01–T28 与正文相互索引；经验材料 E1–E11 汇集全部分析分支及其数据来源。图的矢量源文件、源数据和复现代码随文提供。'
            heads=['A. 理论与统计证明','B. 真实重复测量与校准','C. 与相关研究的关系','D. 参考文献']
        else:
            title='Supplementary Information: Identifying direction and temporal weighting in predictive cognitive maps'
            lead='This file provides the general proofs supporting the main text, complete methods and results for real repeated measurements and known-transformation calibration, and related-work context. Results T01–T28 are cross-referenced in the article; empirical sections E1–E11 collect all analysis branches and their sources. Editable vector figures, source data and reproduction code accompany the paper.'
            heads=['A. Theory and statistical proofs','B. Real repeated measurements and calibration','C. Relation to prior work','D. References']
        front=(f'---\ntitle: "{title}"\nlang: '+('zh-CN' if lang=='zh' else 'en')+
               '\nauthor: "Yongzhang Ji"\naffiliation: "Faculty of Psychology, Southwest University, Chongqing, China"'
               '\norcid: "0009-0006-4111-2000"\ndate: "2026-11-26"\n---\n\n')
        bodies=[no_title((ROOT/f'supplement/Theory_Appendix_{upper}.md').read_text()),
                no_title((ROOT/f'supplement/Supplementary_Empirical_{upper}.md').read_text()),
                no_title((ROOT/f'supplement/Related_Work_{upper}.md').read_text()),refs]
        supp=front+lead+'\n\n'+'\n\n'.join('# '+h+'\n\n'+b for h,b in zip(heads,bodies))+'\n'
        (M/f'Supplementary_{upper}.md').write_text(supp)
        print(f'{upper}: article {len(s):,} characters; supplement {len(supp):,} characters')
    registry=[dict(figure=i,stem=caps[str(i)]['stem'],caption_zh=caps[str(i)]['zh'],caption_en=caps[str(i)]['en']) for i in range(1,9)]
    (ROOT/'source_data/figure_caption_registry.json').write_text(json.dumps(registry,ensure_ascii=False,indent=2)+'\n')

if __name__=='__main__':main()

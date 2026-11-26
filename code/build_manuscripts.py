"""Build four bilingual PDFs and editable TeX from assembled Markdown.

Requires Pandoc, XeLaTeX, standard TeX packages and Noto fonts.
Figures are exported vector PDFs from the native editable SVGs, never raster.
Figures keep their complete scientific captions; Figure 8 opens Discussion.
"""
from pathlib import Path
import argparse, os, re, subprocess, json, shutil

ROOT=Path(__file__).resolve().parents[1]
BUILD=ROOT/'build';BUILD.mkdir(exist_ok=True)
M=ROOT/'manuscripts'

def run(cmd,*,env=None,cwd=ROOT):
    return subprocess.run(cmd,cwd=cwd,env=env,capture_output=True,text=True)

def tex_environment():
    env=dict(os.environ)
    # Some container images ship a complete TeX tree but a dangling ls-R
    # database link. Explicit fallback paths avoid changing the installation.
    if not run(['kpsewhich','article.cls']).stdout.strip():
        dist=Path('/usr/share/texlive/texmf-dist')
        if (dist/'tex/latex/base/article.cls').exists():
            env['TEXINPUTS']='.:'+str(dist/'tex')+'//:'+env.get('TEXINPUTS','')
            env['TFMFONTS']=str(dist/'fonts/tfm')+'//:'+env.get('TFMFONTS','')
            env['VFFONTS']=str(dist/'fonts/vf')+'//:'+env.get('VFFONTS','')
            env['T1FONTS']=str(dist/'fonts/type1')+'//:'+env.get('T1FONTS','')
            env['OPENTYPEFONTS']='/usr/share/texmf/fonts/opentype//:'+str(dist/'fonts/opentype')+'//:'+env.get('OPENTYPEFONTS','')
    found=run(['kpsewhich','xelatex.fmt']).stdout.strip()
    if not found and not (BUILD/'xelatex.fmt').exists():
        old=ROOT/'build/xelatex.fmt.cached'
        if old.exists():shutil.copy2(old,BUILD/'xelatex.fmt')
        else:
            r=run(['xetex','-ini','-etex','-jobname=xelatex','-progname=xelatex',
                   '-interaction=nonstopmode','-halt-on-error','xelatex.ini'],cwd=BUILD,env=env)
            (BUILD/'format_build.log').write_text(r.stdout+r.stderr)
            if r.returncode:raise RuntimeError('XeLaTeX format creation failed')
    env['TEXFORMATS']=str(BUILD)+os.pathsep+env.get('TEXFORMATS','')
    return env

def math_in_prose(s):
    chunks=re.split(r'(\$\$.*?\$\$|\$(?:[^$]|\\\$)*?\$|\\\(.*?\\\)|\\\[.*?\\\])',s,flags=re.S)
    single={'γ':r'\gamma','η':r'\eta','κ':r'\kappa','λ':r'\lambda',
            'Φ':r'\Phi','ω':r'\omega','π':r'\pi','∎':r'\square',
            '⟨':r'\langle','⟩':r'\rangle','‖':r'\Vert','−':'-',
            '⊤':r'\top','α':r'\alpha','σ':r'\sigma','ρ':r'\rho',
            'δ':r'\delta','θ':r'\theta','ε':r'\epsilon','φ':r'\phi',
            '≥':r'\ge','≤':r'\le','≠':r'\ne','→':r'\to','×':r'\times',
            '∞':r'\infty','∈':r'\in','Σ':r'\Sigma','Π':r'\Pi'}
    trans={'ᵀ':r'^\top','⊤':r'^\top','²':'^2','³':'^3','⁴':'^4','⁷':'^7',
           '₀':'_0','₁':'_1','₂':'_2','₃':'_3','₄':'_4'}
    for i in range(0,len(chunks),2):
        x=chunks[i]
        # Keep a Unicode minus and its following number in one math span.
        # Splitting them makes Pandoc interpret a delimiter before a digit as text.
        x=re.sub(r'−(\d+(?:\.\d+)?|\.\d+)',lambda m:r'\('+ '-'+m.group(1)+r'\)',x)
        supers={'⁻':'-','⁺':'+','⁰':'0','¹':'1','²':'2','³':'3','⁴':'4','⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9'}
        x=re.sub(r'(\d+(?:\.\d+)?)([⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹]+)',lambda m:r'\('+m.group(1)+'^{'+''.join(supers[c] for c in m.group(2))+r'}\)',x)
        x=re.sub(r'([A-Za-z])([ᵀ⊤²³⁴⁷₀₁₂₃₄]+)',lambda m:r'\('+m.group(1)+''.join(trans[c] for c in m.group(2))+r'\)',x)
        for a,b in single.items():x=x.replace(a,r'\('+b+r'\)')
        chunks[i]=x
    return ''.join(chunks)

def source_for_pdf(s,lang,supplement):
    if not supplement:
        # Shift body headings explicitly. Pandoc's negative heading shift can
        # otherwise replace YAML title metadata with the first H1 ('Abstract').
        s=re.sub(r'(?m)^(#{1,6}) ',lambda m:'#'*max(1,len(m.group(1))-1)+' ',s)
    def figure(m):
        content=m.group(2).strip()
        im=re.match(r'!\[[^\]]+\]\(\.\./figures/svg/([^\)]+)\.svg\)\s*',content)
        if not im:raise ValueError('Missing figure image')
        stem=im.group(1)
        caption=content[im.end():].strip()
        # A float keeps each vector panel and its full caption together while
        # allowing surrounding prose to fill the preceding page naturally.
        start=('```{=latex}\n\\begin{figure}[!htbp]\n\\centering\n'
               '\\includegraphics[width=\\linewidth,height=0.64\\textheight,keepaspectratio]{figures/pdf/'+stem+'.pdf}\n'
               '\\vspace{3pt}\n\\begin{minipage}{\\linewidth}\n'
               '\\fontsize{9.2}{11.7}\\selectfont\n\\setlength{\\parskip}{3pt}\n```\n\n')
        return start+caption+'\n\n```{=latex}\n\\end{minipage}\n\\end{figure}\n```\n'
    s=re.sub(r'<!--FIGURE_START:(\d+)-->\s*(.*?)\s*<!--FIGURE_END:\1-->',figure,s,flags=re.S)
    if not supplement:
        s=re.sub(r'(?m)^# (Discussion|讨论|Methods|方法)$',lambda m:'\\FloatBarrier\n\n'+m.group(0),s)
        s=re.sub(r'(?m)^# (References|参考文献)$',
                 lambda m:'\\FloatBarrier\n\n\\Needspace{8\\baselineskip}\n\n'+m.group(0)+'\n\n\\fontsize{9}{11}\\selectfont\n\\setlength{\\parskip}{2pt}\n\\setlist[enumerate]{itemsep=2pt}',s)
        s=re.sub(r'(?m)^# (Methods|方法)$',
                 lambda m:'\\Needspace{8\\baselineskip}\n\n'+m.group(0),s)
    if supplement:
        # Each supplement part has a separate entry and page boundary.
        s=re.sub(r'(?m)^# A\. ',lambda m:'\\begingroup\n\\fontsize{8.8}{10.3}\\selectfont\n\\linespread{1.0}\\selectfont\n\\setlength{\\parskip}{0pt}\n\\tableofcontents\n\\endgroup\n\n'+m.group(0),s,count=1)
        s=re.sub(r'(?m)^# ([ABCD]\. .*)$',r'\\clearpage\n\n# \1',s)
        s=s.replace('# C. 已有贡献与本文增量','# C. 已有贡献与本文增量\n\n\\small')
        s=s.replace('# C. Prior contributions and the present increment','# C. Prior contributions and the present increment\n\n\\small')
        def keep_short_table(m):
            block=m.group(0)
            rows=block.strip().splitlines()
            prefix='\\begingroup\n\\small\n\n'
            if len(rows)<=8:
                prefix='\\Needspace{'+str(len(rows)+3)+'\\baselineskip}\n\n'+prefix
            return prefix+block+'\n\n\\endgroup\n'
        s=re.sub(r'(?m)(?:^\|[^\n]+\|\s*\n){3,}',keep_short_table,s)
    if supplement:
        s=re.sub(r'(?m)^# D\. (References|参考文献)$',lambda m:m.group(0)+'\n\n\\fontsize{9}{11}\\selectfont\n\\setlength{\\parskip}{2pt}\n\\setlist[enumerate]{itemsep=2pt}',s)
    # Long machine paths should wrap rather than run into a margin.
    s=re.sub(r'`([^`\n]*[/_][^`\n]*)`',lambda m:r'\path{'+m.group(1)+'}',s)
    # Keep a display with the paragraph that introduces it. This also prevents
    # a floating figure from separating an introduction from its equation.
    def joined_display(m):
        return ('```{=latex}\n\\par\n\\noindent\\begin{minipage}{\\linewidth}\n```\n\n'+m.group(1)+'\n\n'+m.group(2)+
                '\n\n```{=latex}\n\\end{minipage}\\par\n```')
    s=re.sub(r'(?m)^((?![#|$`\\])\S[^\n]*)\n\n(\$\$\n.*?\n\$\$)',joined_display,s,flags=re.S)
    return math_in_prose(s)

def build(name,env):
    lang='ZH' if name.endswith('_ZH') else 'EN'
    supplement=name.startswith('Supplementary')
    s=(M/(name+'.md')).read_text()
    temp=BUILD/(name+'.render.md');temp.write_text(source_for_pdf(s,lang,supplement))
    tex=M/(name+'.tex')
    kind=('补充材料' if supplement else '预印本') if lang=='ZH' else ('Supplementary Information' if supplement else 'Preprint')
    running='认知预测地图的机制识别' if lang=='ZH' else 'Identifying mechanisms in predictive maps'
    cmd=['pandoc',str(temp),'--from=markdown+tex_math_dollars+tex_math_single_backslash+raw_tex+pipe_tables+implicit_figures+autolink_bare_uris',
         '--to=latex','--standalone','--template='+str(ROOT/'code/paper_template.tex'),
         '-V','document_kind='+kind,'-V','running_title='+running,'--resource-path='+str(ROOT),'-o',str(tex)]
    if lang=='ZH':cmd+=['-V','chinese=true']
    r=run(cmd);(BUILD/(name+'.pandoc.log')).write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Pandoc failed: '+r.stderr[-1500:])
    # Use real Unicode punctuation, independently of font-specific TeX
    # ligatures. Preserve math, URLs, paths and verbatim code unchanged.
    ts=tex.read_text()
    ts=re.sub(r'(?m)^(\\subsection\{)',r'\\Needspace{7\\baselineskip}\n\1',ts)
    ts=re.sub(r'(?m)^(\\subsubsection\{)',r'\\Needspace{6\\baselineskip}\n\1',ts)
    spans=re.split(r'(\\\(.*?\\\)|\\\[.*?\\\]|\\(?:path|url|href)\{[^}]*\}|\\begin\{verbatim\}.*?\\end\{verbatim\})',ts,flags=re.S)
    for i in range(0,len(spans),2):
        spans[i]=spans[i].replace('---','—').replace('--','–').replace('``','“').replace("''",'”')
    tex.write_text(''.join(spans))
    previous_aux=None
    for k in range(4):
        r=run(['xelatex','-interaction=nonstopmode','-halt-on-error','-file-line-error',
               '-output-directory='+str(BUILD),str(tex)],env=env)
        (BUILD/(name+f'.pass{k+1}.log')).write_text(r.stdout+r.stderr)
        if r.returncode:raise RuntimeError(f'XeLaTeX failed for {name}: '+(r.stdout+r.stderr)[-2500:])
        current_aux=''.join((BUILD/(name+ext)).read_text() if (BUILD/(name+ext)).exists() else '' for ext in ['.aux','.toc'])
        if k>=2 and current_aux==previous_aux:break
        previous_aux=current_aux
    shutil.copy2(BUILD/(name+'.pdf'),M/(name+'.pdf'))
    log=(BUILD/(name+'.log')).read_text(errors='replace')
    report={'name':name,'bytes':(M/(name+'.pdf')).stat().st_size,'overfull_boxes':log.count('Overfull'),
            'missing_glyph_lines':[l for l in log.splitlines() if 'Missing character:' in l],
            'undefined_references':'undefined references' in log}
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('names',nargs='*');args=ap.parse_args()
    names=args.names or ['Manuscript_ZH','Manuscript_EN','Supplementary_ZH','Supplementary_EN']
    env=tex_environment();results=[]
    for n in names:
        report=build(n,env);results.append(report);print(json.dumps(report,ensure_ascii=False),flush=True)
    report_file=BUILD/'build_report.json'
    previous=json.loads(report_file.read_text()) if report_file.exists() else []
    merged={r['name']:r for r in previous}
    merged.update({r['name']:r for r in results})
    report_file.write_text(json.dumps(list(merged.values()),ensure_ascii=False,indent=2))

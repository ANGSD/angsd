"""Turn a GitHub wiki checkout into an MkDocs Material site."""
import os,re,shutil,sys,yaml
WIKI,OUT=sys.argv[1],sys.argv[2]
WIKI_URL=os.environ.get('WIKI_URL','https://github.com/ANGSD/angsd/wiki')
REPO_URL=os.environ.get('REPO_URL','https://github.com/ANGSD/angsd')
DOCS=os.path.join(OUT,'docs')
shutil.rmtree(DOCS,ignore_errors=True); os.makedirs(DOCS)

pages={}   # slug -> filename in docs/
for f in sorted(os.listdir(WIKI)):
    if f.endswith('.md') and not f.startswith('_'):
        slug=f[:-3].replace(' ','-')
        pages[slug]='index.md' if slug=='Home' else slug+'.md'
    elif f.lower().endswith('.png'):
        shutil.copy(os.path.join(WIKI,f),os.path.join(DOCS,f))

def fixlinks(txt):
    def rep(m):
        label,tgt=m.group(1),m.group(2)
        if re.match(r'^(https?:|mailto:|#|/)',tgt): return m.group(0)
        base,_,anchor=tgt.partition('#')
        if base in pages:
            return '[%s](%s%s)'%(label,pages[base],'#'+anchor if anchor else '')
        return m.group(0)
    return re.sub(r'\[([^\]]*)\]\(([^)\s]+)\)',rep,txt)

for f in sorted(os.listdir(WIKI)):
    if not f.endswith('.md') or f.startswith('_'): continue
    slug=f[:-3].replace(' ','-')
    txt=fixlinks(open(os.path.join(WIKI,f)).read())
    if not re.match(r'^\s*#\s',txt):                     # give every page an H1
        txt='# %s\n\n'%f[:-3]+txt
    link=WIKI_URL if slug=='Home' else WIKI_URL+'/'+slug
    txt+=('\n\n---\n\n*This page is a snapshot. The live, always-current version is on the '
          '[ANGSD wiki](%s) \u2014 edit it there and this site rebuilds automatically.*\n'%link)
    open(os.path.join(DOCS,pages[slug]),'w').write(txt)

# nav from the wiki's own _Sidebar.md
nav=[{'ANGSD':'index.md'}]; sec=None; items=[]
for ln in open(os.path.join(WIKI,'_Sidebar.md')):
    ln=ln.rstrip()
    m=re.match(r'^\*\*(.+?)\*\*\s*$',ln)
    if m:
        if sec and items: nav.append({sec:items})
        sec=re.sub(r'\[|\]\(.*?\)','',m.group(1)); items=[]; continue
    m=re.match(r'^- \[(.+?)\]\(([^)#]+)(#[^)]*)?\)',ln)
    if m and sec:
        label,tgt=m.group(1),m.group(2)
        if tgt=='Home': tgt='index'
        f=pages.get(tgt)
        # nav can't point at an anchor inside a file; use the page's URL instead
        if f and m.group(3): items.append({label:('' if f=='index.md' else f[:-3]+'/')+m.group(3)})
        elif f: items.append({label:f})
if sec and items: nav.append({sec:items})
listed={list(i.values())[0].split('#')[0] for s in nav for i in (list(s.values())[0] if isinstance(list(s.values())[0],list) else [s])}
rest=sorted(f for f in pages.values() if f not in listed and f!='index.md')
if rest: nav.append({'Other pages':rest})

cfg={
 'site_name':'ANGSD',
 'site_description':'Analysis of Next Generation Sequencing Data',
 'repo_url':REPO_URL,
 'edit_uri_template':'wiki/{path_noext}',
 'theme':{'name':'material','features':['navigation.sections','navigation.top','content.code.copy','search.suggest','toc.follow'],
   'palette':[{'media':'(prefers-color-scheme: light)','scheme':'default','primary':'teal','toggle':{'icon':'material/weather-night','name':'Dark mode'}},
              {'media':'(prefers-color-scheme: dark)','scheme':'slate','primary':'teal','toggle':{'icon':'material/weather-sunny','name':'Light mode'}}]},
 'markdown_extensions':['tables','admonition','attr_list','md_in_html','toc',
   {'pymdownx.arithmatex':{'generic':True}},
   {'pymdownx.superfences':{'custom_fences':[{'name':'mermaid','class':'mermaid','format':'!!python/name:pymdownx.superfences.fence_code_format'}]}},
   'pymdownx.details','pymdownx.highlight','pymdownx.inlinehilite'],
 'extra_javascript':['js/mathjax.js','https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-svg.js'],
 'nav':nav,
 'validation':{'nav':{'not_found':'info'}},   # the anchor URLs above aren't files
}
os.makedirs(os.path.join(DOCS,'js'),exist_ok=True)
open(os.path.join(DOCS,'js','mathjax.js'),'w').write(
 'window.MathJax={tex:{inlineMath:[["\\\\(","\\\\)"]],'
 'displayMath:[["\\\\[","\\\\]"]]},svg:{fontCache:"local"},'
 'options:{skipHtmlTags:["script","noscript","style","textarea","pre","code"]}};')
txt=yaml.dump(cfg,sort_keys=False,default_flow_style=False,allow_unicode=True)
txt=txt.replace("'!!python/name:pymdownx.superfences.fence_code_format'","!!python/name:pymdownx.superfences.fence_code_format")
open(os.path.join(OUT,'mkdocs.yml'),'w').write(txt)
print("docs/: %d pages, %d images"%(len(pages),len([f for f in os.listdir(DOCS) if f.endswith('.png')])))
print("nav sections: %d"%len(nav))

#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import fitz

def place(page: fitz.Page, path: Path, box):
    src = fitz.open(path)
    page.show_pdf_page(fitz.Rect(*box), src, 0, keep_proportion=True)
    src.close()

def main():
    p=argparse.ArgumentParser(); p.add_argument("root",type=Path); a=p.parse_args()
    f=a.root.resolve()/"results"/"figures"
    doc=fitz.open(); page=doc.new_page(width=1728,height=1220)
    panels=[
        ("Figure5A_GDSC2_volcano.pdf",(0,0,560,380)),
        ("Figure5B_CTRPv2_volcano.pdf",(560,0,1120,380)),
        ("Figure5C_PRISM_volcano.pdf",(1120,0,1728,380)),
        ("Figure5D_significant_pair_membership.pdf",(0,380,510,780)),
        ("Figure5E_replicated_association_heatmap.pdf",(510,375,1160,820)),
        ("Figure5F_TF_drug_integration.pdf",(1160,380,1728,820)),
        ("Figure5G_PDX_validation.pdf",(230,820,865,1220)),
        ("Figure5H_PDX_waterfall.pdf",(865,820,1500,1220)),
    ]
    for name,box in panels: place(page,f/name,box)
    out=f/"Figure5_complete.pdf"; doc.save(out,garbage=4,deflate=True); doc.close()
    for names,outname,size in [
        (["SupplementaryFigure8A_cellline_mapping.pdf","SupplementaryFigure8B_drug_coverage.pdf","SupplementaryFigure8C_classifier_validation.pdf"],"SupplementaryFigure8_complete.pdf",(1400,620)),
        (["SupplementaryFigure9A_PDX_drug_coverage.pdf","SupplementaryFigure9B_PDX_association_screen.pdf"],"SupplementaryFigure9_complete.pdf",(1200,620)),
    ]:
        d=fitz.open(); pg=d.new_page(width=size[0],height=size[1]); n=len(names)
        for i,name in enumerate(names): place(pg,f/name,(i*size[0]/n,0,(i+1)*size[0]/n,size[1]))
        d.save(f/outname,garbage=4,deflate=True); d.close()
    for name in ("Figure5_complete.pdf","SupplementaryFigure8_complete.pdf","SupplementaryFigure9_complete.pdf"):
        d=fitz.open(f/name); pix=d[0].get_pixmap(matrix=fitz.Matrix(1.4,1.4),alpha=False); pix.save((f/name).with_suffix(".png")); d.close()
    return 0
if __name__=="__main__": raise SystemExit(main())


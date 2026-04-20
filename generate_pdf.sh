#!/bin/bash

# pandocがインストールされているか確認
if ! command -v pandoc &> /dev/null; then
    echo "pandocがインストールされていません"
    echo "macOSの場合: brew install pandoc"
    echo "または: brew install basictex pandoc"
    exit 1
fi

# PDF生成
pandoc dilithium2_report_ja.md \
    -o dilithium2_report_ja.pdf \
    --pdf-engine=xelatex \
    -V CJKmainfont="Hiragino Mincho ProN" \
    -V geometry:margin=1in \
    --toc \
    --toc-depth=2 \
    --number-sections \
    -V linkcolor:blue \
    -V urlcolor:blue

if [ $? -eq 0 ]; then
    echo "✅ PDFが正常に生成されました: dilithium2_report_ja.pdf"
    ls -lh dilithium2_report_ja.pdf
else
    echo "❌ PDF生成に失敗しました"
fi

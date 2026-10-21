"""Export the generated CSV as an editable landscape Word table (stdlib only)."""
from pathlib import Path
import csv
from html import escape
from zipfile import ZipFile, ZIP_DEFLATED
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
with (HERE / 'ADNI_paired_omics_characteristics.csv').open(encoding='utf-8') as handle:
    rows = list(csv.reader(handle))
notes = (HERE / 'ADNI_paired_omics_characteristics.md').read_text(encoding='utf-8').splitlines()[-1]
assert [row[0] for row in rows].count('APOE ε4 status, n (%)') == 1
assert '<U+' not in str(rows) + notes

def paragraph(text, bold=False, size=18):
    return ('<w:p><w:pPr><w:spacing w:after="40"/></w:pPr><w:r><w:rPr>'
            '<w:rFonts w:ascii="Arial" w:hAnsi="Arial"/>'
            + ('<w:b/>' if bold else '') + f'<w:sz w:val="{size}"/>'
            + '</w:rPr><w:t xml:space="preserve">' + escape(text) + '</w:t></w:r></w:p>')

widths = [2500] + [1250] * (len(rows[0]) - 3) + [1400, 1400]
table_rows = []
for i, row in enumerate(rows):
    is_section = i > 0 and not any(row[1:-2])
    cells = []
    for width, text in zip(widths, row):
        shade = '<w:shd w:fill="EEF2F5"/>' if i == 0 or is_section else ''
        cells.append(f'<w:tc><w:tcPr><w:tcW w:w="{width}" w:type="dxa"/>{shade}</w:tcPr>'
                     + paragraph(text, bold=i == 0 or is_section) + '</w:tc>')
    table_rows.append('<w:tr>' + ('<w:trPr><w:tblHeader/></w:trPr>' if i == 0 else '')
                      + ''.join(cells) + '</w:tr>')
table = (f'<w:tbl><w:tblPr><w:tblW w:w="{sum(widths)}" w:type="dxa"/>'
         '<w:tblBorders><w:top w:val="single" w:sz="6"/>'
         '<w:bottom w:val="single" w:sz="6"/>'
         '<w:insideH w:val="single" w:sz="3" w:color="DDDDDD"/></w:tblBorders>'
         '</w:tblPr><w:tblGrid>'
         + ''.join(f'<w:gridCol w:w="{w}"/>' for w in widths)
         + '</w:tblGrid>' + ''.join(table_rows) + '</w:tbl>')
document = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>'
            + paragraph('Table 1. Baseline characteristics of the ADNI paired-omics cohort by diagnosis and MCI subtype', True, 22)
            + table + paragraph(notes, size=17)
            + '<w:sectPr><w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/>'
            '<w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720"/>'
            '</w:sectPr></w:body></w:document>')
ET.fromstring(document)
content_types = ('<?xml version="1.0" encoding="UTF-8"?>'
                 '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                 '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                 '<Default Extension="xml" ContentType="application/xml"/>'
                 '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
                 '</Types>')
relationships = ('<?xml version="1.0" encoding="UTF-8"?>'
                 '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                 '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>'
                 '</Relationships>')
output = HERE / 'ADNI_paired_omics_characteristics.docx'
with ZipFile(output, 'w', ZIP_DEFLATED) as archive:
    archive.writestr('[Content_Types].xml', content_types)
    archive.writestr('_rels/.rels', relationships)
    archive.writestr('word/document.xml', document)
with ZipFile(output) as archive:
    assert archive.testzip() is None
print(output)

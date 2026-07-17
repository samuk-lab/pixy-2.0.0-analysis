#!/usr/bin/env bash
# figure 1 — hand-built pixy capabilities + development schematic (run under WSL)
# writes Figure_1_pixy_overview.svg (vector master) and figs/Figure_pixy_overview.pdf
# svg -> pdf via rsvg-convert / inkscape / cairosvg, else Windows Edge over WSL interop
#   bash Figure_1_pixy_overview.sh
# restyle: edit the colors block below; relayout: edit the svg heredoc

set -euo pipefail

##########
# colors — edit to restyle
##########
BG="#ffffff"            # page background
TEXT_PRIMARY="#2C2C2A"  # column headers / neutral titles
TEXT_SECONDARY="#5F5E5A" # neutral captions
ARROW="#888780"         # connector arrows
DIVIDER="#B4B2A9"        # horizontal rule between the two bands

# per category: box fill / box stroke / title text / subtitle text
# inputs band
TEAL_FILL="#E1F5EE";   TEAL_STROKE="#0F6E56";   TEAL_TITLE="#085041";   TEAL_SUB="#0F6E56"
# pixy core band
PURPLE_FILL="#EEEDFE"; PURPLE_STROKE="#534AB7"; PURPLE_TITLE="#3C3489"; PURPLE_SUB="#534AB7"
# outputs (per-statistic) band
CORAL_FILL="#FAECE7";  CORAL_STROKE="#993C1D";  CORAL_TITLE="#712B13";  CORAL_SUB="#993C1D"
# development & distribution band
BLUE_FILL="#E6F1FB";   BLUE_STROKE="#185FA5";   BLUE_TITLE="#0C447C";   BLUE_SUB="#185FA5"
# neutral note box (aggregation hint)
GRAY_FILL="#F1EFE8";   GRAY_STROKE="#5F5E5A";   GRAY_TEXT="#444441"
# "new" badge
AMBER_FILL="#FAC775";  AMBER_STROKE="#854F0B";  AMBER_TEXT="#633806"

# font stack for all text
FONT="'Helvetica Neue',Arial,sans-serif"

##########
# paths
##########
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SVG="$SCRIPT_DIR/Figure_1_pixy_overview.svg"
PDF_DIR="$SCRIPT_DIR/figs"
PDF="$PDF_DIR/Figure_pixy_overview.pdf"
mkdir -p "$PDF_DIR"

##########
# emit the svg (single source of truth)
##########
# colours from the variables above; special glyphs use XML numeric entities to
# keep the body pure ASCII
cat > "$SVG" <<EOF
<svg width="680" height="420" viewBox="0 0 680 420" xmlns="http://www.w3.org/2000/svg" role="img">
<title>pixy capabilities and development flowchart</title>
<desc>A left-to-right pixy analysis pipeline (inputs, pixy core, outputs) plus a parallel development and distribution band. Capabilities added since the original pixy release are badged "new".</desc>
<style>
  text{font-family:${FONT}}
  .t{font-weight:400;font-size:14px;fill:${TEXT_PRIMARY}}
  .th{font-weight:500;font-size:14px;fill:${TEXT_PRIMARY}}
  .ts{font-weight:400;font-size:12px;fill:${TEXT_SECONDARY}}
  .c-teal>rect{fill:${TEAL_FILL};stroke:${TEAL_STROKE};stroke-width:.5}
  .c-teal .th{fill:${TEAL_TITLE}}.c-teal .ts{fill:${TEAL_SUB}}
  .c-purple>rect{fill:${PURPLE_FILL};stroke:${PURPLE_STROKE};stroke-width:.5}
  .c-purple .th{fill:${PURPLE_TITLE}}.c-purple .ts{fill:${PURPLE_SUB}}
  .c-coral>rect{fill:${CORAL_FILL};stroke:${CORAL_STROKE};stroke-width:.5}
  .c-coral .th{fill:${CORAL_TITLE}}.c-coral .ts{fill:${CORAL_SUB}}
  .c-blue>rect{fill:${BLUE_FILL};stroke:${BLUE_STROKE};stroke-width:.5}
  .c-blue .th{fill:${BLUE_TITLE}}.c-blue .ts{fill:${BLUE_SUB}}
  .c-gray>rect{fill:${GRAY_FILL};stroke:${GRAY_STROKE};stroke-width:.5}
  .c-gray .ts{fill:${GRAY_TEXT}}
  .c-amber>rect{fill:${AMBER_FILL};stroke:${AMBER_STROKE};stroke-width:.5}
  .c-amber .ts{fill:${AMBER_TEXT}}
</style>
<defs><marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M2 1L8 5L2 9" fill="none" stroke="${ARROW}" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/></marker></defs>

<rect x="0" y="0" width="680" height="420" fill="${BG}"/>

<text class="th" x="130" y="34" text-anchor="middle">Inputs</text>
<text class="th" x="385" y="34" text-anchor="middle">pixy core</text>
<text class="th" x="590" y="34" text-anchor="middle">Outputs</text>

<g class="c-teal"><rect x="24" y="56" width="212" height="104" rx="6"/>
<text class="th" x="36" y="78">Genotype VCF</text>
<text class="ts" x="36" y="94">bgzip + tabix/CSI index</text>
<text class="ts" x="36" y="116">Callable-site source:</text>
<text class="ts" x="44" y="134">&#8226; all-sites VCF</text>
<text class="ts" x="44" y="152">&#8226; GVCF blocks</text>
</g>
<g class="c-amber"><rect x="186" y="84" width="34" height="15" rx="7"/><text class="ts" x="203" y="95" text-anchor="middle">new</text></g>

<g class="c-teal"><rect x="24" y="172" width="212" height="44" rx="6"/>
<text class="th" x="36" y="194">Populations file</text>
<text class="ts" x="36" y="210">sample &#8594; population (TSV)</text>
</g>

<g class="c-teal"><rect x="24" y="228" width="212" height="62" rx="6"/>
<text class="th" x="36" y="250">Optional companions</text>
<text class="ts" x="44" y="268">BED &#8594; custom windows</text>
<text class="ts" x="44" y="284">sites &#8594; site subset</text>
</g>
<g class="c-amber"><rect x="196" y="258" width="34" height="15" rx="7"/><text class="ts" x="213" y="269" text-anchor="middle">new</text></g>
<g class="c-amber"><rect x="196" y="274" width="34" height="15" rx="7"/><text class="ts" x="213" y="285" text-anchor="middle">new</text></g>

<line x1="240" y1="174" x2="296" y2="174" stroke="${ARROW}" stroke-width="1.5" marker-end="url(#arrow)"/>

<g class="c-purple"><rect x="300" y="56" width="170" height="60" rx="6"/>
<text class="th" x="312" y="78">Windowing</text>
<text class="ts" x="312" y="94">window / BED / region</text>
<text class="ts" x="312" y="110">chromosome subset</text>
</g>
<g class="c-purple"><rect x="300" y="126" width="170" height="44" rx="6"/>
<text class="th" x="312" y="148">Statistics</text>
<text class="ts" x="312" y="164">select via --stats</text>
</g>
<g class="c-purple"><rect x="300" y="180" width="170" height="112" rx="6"/>
<text class="th" x="312" y="202">Options</text>
<text class="ts" x="312" y="224">Hudson FST</text>
<text class="ts" x="312" y="244">multiallelic SNPs</text>
<text class="ts" x="312" y="264">ploidy handling</text>
<text class="ts" x="312" y="284">multicore</text>
</g>
<g class="c-amber"><rect x="430" y="214" width="34" height="15" rx="7"/><text class="ts" x="447" y="225" text-anchor="middle">new</text></g>
<g class="c-amber"><rect x="430" y="234" width="34" height="15" rx="7"/><text class="ts" x="447" y="245" text-anchor="middle">new</text></g>
<g class="c-amber"><rect x="430" y="254" width="34" height="15" rx="7"/><text class="ts" x="447" y="265" text-anchor="middle">new</text></g>
<g class="c-amber"><rect x="430" y="274" width="34" height="15" rx="7"/><text class="ts" x="447" y="285" text-anchor="middle">new</text></g>

<line x1="474" y1="174" x2="510" y2="174" stroke="${ARROW}" stroke-width="1.5" marker-end="url(#arrow)"/>

<text class="ts" x="516" y="48">one TSV per statistic</text>
<g class="c-coral"><rect x="516" y="56" width="148" height="32" rx="5"/>
<text class="th" x="528" y="76">&#960;</text><text class="ts" x="546" y="76">within-pop</text>
</g>
<g class="c-coral"><rect x="516" y="92" width="148" height="32" rx="5"/>
<text class="th" x="528" y="112">dxy</text><text class="ts" x="554" y="112">between-pop</text>
</g>
<g class="c-coral"><rect x="516" y="128" width="148" height="32" rx="5"/>
<text class="th" x="528" y="148">FST</text><text class="ts" x="556" y="148">WC / Hudson</text>
</g>
<g class="c-coral"><rect x="516" y="164" width="148" height="32" rx="5"/>
<text class="th" x="528" y="184">Watterson &#952;</text>
</g>
<g class="c-amber"><rect x="620" y="172" width="34" height="15" rx="7"/><text class="ts" x="637" y="183" text-anchor="middle">new</text></g>
<g class="c-coral"><rect x="516" y="200" width="148" height="32" rx="5"/>
<text class="th" x="528" y="220">Tajima's D</text>
</g>
<g class="c-amber"><rect x="620" y="208" width="34" height="15" rx="7"/><text class="ts" x="637" y="219" text-anchor="middle">new</text></g>

<g class="c-gray"><rect x="516" y="238" width="148" height="52" rx="5"/>
<text class="ts" x="528" y="256">raw count / component</text>
<text class="ts" x="528" y="271">columns &#8594; exact</text>
<text class="ts" x="528" y="286">window aggregation</text>
</g>

<line x1="24" y1="308" x2="664" y2="308" stroke="${DIVIDER}" stroke-width="0.5"/>
<text class="th" x="24" y="328">Development &amp; distribution</text>

<g class="c-blue"><rect x="24" y="338" width="150" height="58" rx="6"/>
<text class="th" x="36" y="360">Contribute</text>
<text class="ts" x="36" y="378">GitHub fork + PR</text>
<text class="ts" x="36" y="392">poetry install</text>
</g>
<g class="c-amber"><rect x="134" y="345" width="34" height="15" rx="7"/><text class="ts" x="151" y="356" text-anchor="middle">new</text></g>
<line x1="174" y1="367" x2="198" y2="367" stroke="${ARROW}" stroke-width="1.5" marker-end="url(#arrow)"/>
<g class="c-blue"><rect x="198" y="338" width="150" height="58" rx="6"/>
<text class="th" x="210" y="360">Local checks</text>
<text class="ts" x="210" y="378">ruff &#183; mypy &#183; pytest</text>
</g>
<g class="c-amber"><rect x="308" y="345" width="34" height="15" rx="7"/><text class="ts" x="325" y="356" text-anchor="middle">new</text></g>
<line x1="348" y1="367" x2="372" y2="367" stroke="${ARROW}" stroke-width="1.5" marker-end="url(#arrow)"/>
<g class="c-blue"><rect x="372" y="338" width="150" height="58" rx="6"/>
<text class="th" x="384" y="360">CI on every PR</text>
<text class="ts" x="384" y="378">Python 3.10&#8211;3.14</text>
</g>
<g class="c-amber"><rect x="482" y="345" width="34" height="15" rx="7"/><text class="ts" x="499" y="356" text-anchor="middle">new</text></g>
<line x1="522" y1="367" x2="546" y2="367" stroke="${ARROW}" stroke-width="1.5" marker-end="url(#arrow)"/>
<g class="c-blue"><rect x="546" y="338" width="118" height="58" rx="6"/>
<text class="th" x="558" y="360">Distribute</text>
<text class="ts" x="558" y="378">conda-forge</text>
</g>
<g class="c-amber"><rect x="624" y="345" width="34" height="15" rx="7"/><text class="ts" x="641" y="356" text-anchor="middle">new</text></g>

<g class="c-amber"><rect x="24" y="404" width="34" height="15" rx="7"/><text class="ts" x="41" y="415" text-anchor="middle">new</text></g>
<text class="ts" x="66" y="415">added since the original pixy release (Korunes &amp; Samuk 2021)</text>
</svg>
EOF
echo "wrote $SVG"

##########
# convert svg -> pdf: native tool if present, else Windows Edge via interop
##########
if command -v rsvg-convert >/dev/null 2>&1; then
  rsvg-convert -f pdf -o "$PDF" "$SVG"
  echo "wrote $PDF (rsvg-convert)"
elif command -v inkscape >/dev/null 2>&1; then
  inkscape "$SVG" --export-type=pdf --export-filename="$PDF" >/dev/null 2>&1
  echo "wrote $PDF (inkscape)"
elif command -v cairosvg >/dev/null 2>&1; then
  cairosvg "$SVG" -o "$PDF"
  echo "wrote $PDF (cairosvg)"
else
  EDGE="/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
  [ -f "$EDGE" ] || EDGE="/mnt/c/Program Files/Microsoft/Edge/Application/msedge.exe"
  if [ ! -f "$EDGE" ]; then
    echo "ERROR: no SVG->PDF converter found (rsvg-convert/inkscape/cairosvg) and Microsoft Edge is not available for fallback." >&2
    echo "Install one, e.g.:  sudo apt-get install -y librsvg2-bin" >&2
    exit 1
  fi
  # wrap svg in a print-sized html page so Edge prints 1:1, no margins
  HTML="$(mktemp --suffix=.html)"
  {
    printf '%s' '<!DOCTYPE html><html><head><meta charset="utf-8"><style>@page{size:680px 420px;margin:0}html,body{margin:0;padding:0}svg{display:block}</style></head><body>'
    cat "$SVG"
    printf '%s' '</body></html>'
  } > "$HTML"
  HTML_WIN="$(wslpath -m "$HTML")"
  PDF_WIN="$(wslpath -m "$PDF")"
  "$EDGE" --headless=new --disable-gpu --no-pdf-header-footer "--print-to-pdf=${PDF_WIN}" "file:///${HTML_WIN}" >/dev/null 2>&1 || true
  rm -f "$HTML"
  [ -f "$PDF" ] && echo "wrote $PDF (Edge headless via WSL interop)" || { echo "ERROR: PDF was not produced." >&2; exit 1; }
fi

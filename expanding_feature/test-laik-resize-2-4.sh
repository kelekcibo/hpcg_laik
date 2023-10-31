#!/bin/sh
timeout() { perl -e 'alarm shift; exec @ARGV' "$@"; }
timeout 2 ./tcp2run -n 2 -s 2 ./../bin/xhpcg > test-laik-resize-2-4.out
# cmp test-laik-resize-2-4.out "$(dirname -- "${0}")/test-laik-resize-2-4.expected"
echo "Testing expanding the world from size 2 to 4"
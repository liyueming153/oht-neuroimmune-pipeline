#!/bin/bash
# FinnGen H7_GLAUCOMA 看门狗下载: 低速自动断线重连抢初始高速窗口
URL="https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_H7_GLAUCOMA.gz"
OUT=/tmp/finngen_R10_H7_GLAUCOMA.gz
TARGET=805690085
for i in $(seq 1 60); do
  sz=$(stat -c%s "$OUT" 2>/dev/null || echo 0)
  if [ "$sz" -ge "$TARGET" ]; then echo "COMPLETE size=$sz"; break; fi
  echo "attempt $i: have $((sz/1048576))MB, resuming..."
  curl -s -C - --speed-time 25 --speed-limit 20000 "$URL" -o "$OUT"
  sleep 2
done
sz=$(stat -c%s "$OUT" 2>/dev/null || echo 0)
echo "FINAL size=$sz"

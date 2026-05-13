#!/bin/bash
# 监控create_hurdler_lookup.py的进度

echo "================================"
echo "HURDLER Lookup Creation Monitor"
echo "================================"
echo ""

# 检查进程是否在运行
PID=$(ps aux | grep "create_hurdler_lookup.py" | grep -v grep | awk '{print $2}')

if [ -z "$PID" ]; then
    echo "❌ Process not running"
    echo ""
    echo "Checking for output files..."
    ls -lh output/hurdler*.pkl 2>/dev/null || echo "No pkl files found"
    echo ""
    echo "Last 20 lines of log:"
    tail -20 create_lookup_nohup.log 2>/dev/null || tail -20 create_lookup.log 2>/dev/null || echo "No log file found"
else
    echo "✓ Process running (PID: $PID)"
    echo ""
    echo "Latest progress:"
    tail -5 create_lookup_nohup.log 2>/dev/null || tail -5 create_lookup.log 2>/dev/null
    echo ""
    echo "Files created so far:"
    ls -lh output/hurdler*.pkl output/hurdler*.csv 2>/dev/null || echo "No output files yet"
fi

echo ""
echo "================================"

#!/bin/bash
# 监控lookup dictionary创建进度

echo "Monitoring HURDLER Lookup Dictionary Creation"
echo "=============================================="
echo ""

while true; do
    clear
    echo "=== $(date) ==="
    echo ""
    
    # 检查进程
    if ps aux | grep -v grep | grep "create_hurdler_lookup_from_df2.py" > /dev/null; then
        echo "✓ Process is running"
        ps aux | grep -v grep | grep "create_hurdler_lookup_from_df2.py" | awk '{printf "  PID: %s, CPU: %s%%, MEM: %s%%\n", $2, $3, $4}'
    else
        echo "✗ Process not found"
    fi
    
    echo ""
    echo "=== Latest Log Output ==="
    tail -20 create_lookup_from_df2.log 2>/dev/null || echo "No log yet"
    
    echo ""
    echo "=== Output Files ==="
    if [ -f "output/hurdler_lookup_dict.pkl" ]; then
        ls -lh output/hurdler_lookup_dict.pkl | awk '{printf "  hurdler_lookup_dict.pkl: %s\n", $5}'
    fi
    if [ -f "output/hurdler_lookup_lite.pkl" ]; then
        ls -lh output/hurdler_lookup_lite.pkl | awk '{printf "  hurdler_lookup_lite.pkl: %s\n", $5}'
    fi
    
    echo ""
    echo "Press Ctrl+C to stop monitoring"
    
    sleep 30
done

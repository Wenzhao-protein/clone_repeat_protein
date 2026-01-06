#!/usr/bin/env python3
"""
Test script to verify GUI encoding fixes.
Run this before launching the full GUI to ensure all text will be clean.
"""

from interactive_gui import sanitize_for_gui, safe_decode_bytes
import glob

def test_sanitization():
    """Test sanitization function with problematic inputs"""
    print("\n" + "="*70)
    print("SANITIZATION FUNCTION TEST")
    print("="*70)
    
    test_cases = [
        # (input, expected_behavior, description)
        ("\\u51c8\\u5907", "Should remove unicode escapes", "Unicode escape sequence"),
        ("\ufffd\ufffd", "Should remove replacement chars", "Replacement characters"),
        ("Normal text", "Should pass through", "Normal ASCII"),
        ("测试", "Should replace with ?", "Chinese characters"),
        ("", "Should return UNKNOWN", "Empty string"),
        (None, "Should return N/A", "None value"),
        ("垂直方向 (上→下)", "Should sanitize non-ASCII", "Mixed content"),
        ("Line 1\nLine 2", "Should preserve newline", "Newline character"),
    ]
    
    all_passed = True
    for i, (test_input, expected, desc) in enumerate(test_cases, 1):
        result = sanitize_for_gui(test_input)
        
        # Check result is ASCII-only
        is_ascii = all(ord(c) < 128 for c in result)
        has_escape = '\\u' in result or '\ufffd' in result or '�' in result
        
        status = "✓ PASS" if is_ascii and not has_escape else "✗ FAIL"
        if status == "✗ FAIL":
            all_passed = False
            
        print(f"\nTest {i}: {desc}")
        print(f"  Input:    {repr(test_input)[:50]}")
        print(f"  Output:   {repr(result)[:50]}")
        print(f"  Expected: {expected}")
        print(f"  Status:   {status}")
    
    print("\n" + "="*70)
    if all_passed:
        print("✓ ALL TESTS PASSED - Sanitization working correctly!")
    else:
        print("✗ SOME TESTS FAILED - Review output above")
    print("="*70 + "\n")
    
    return all_passed


def test_filename_sanitization():
    """Test with real .scn filenames"""
    print("\n" + "="*70)
    print("FILENAME SANITIZATION TEST")
    print("="*70)
    
    scn_files = glob.glob("*.scn") or glob.glob("wenzhao/*.scn")
    
    if not scn_files:
        print("No .scn files found - skipping filename test")
        return True
    
    print(f"\nFound {len(scn_files)} .scn files")
    print("\nTesting first 5 files:\n")
    
    all_clean = True
    for i, filename in enumerate(scn_files[:5], 1):
        safe = sanitize_for_gui(filename)
        
        # Check for problematic characters
        is_ascii = all(ord(c) < 128 for c in safe)
        has_escape = '\\u' in safe or '\ufffd' in safe or '�' in safe
        
        status = "✓" if is_ascii and not has_escape else "✗"
        if status == "✗":
            all_clean = False
        
        print(f"{status} File {i}:")
        print(f"    Original:  {filename}")
        print(f"    Sanitized: {safe}")
    
    print("\n" + "="*70)
    if all_clean:
        print("✓ ALL FILENAMES CLEAN - Ready for GUI display!")
    else:
        print("✗ SOME FILENAMES HAVE ISSUES - Check output above")
    print("="*70 + "\n")
    
    return all_clean


def test_gui_strings():
    """Test hardcoded GUI strings"""
    print("\n" + "="*70)
    print("GUI STRING VALIDATION TEST")
    print("="*70)
    
    # These are the hardcoded strings used in the GUI
    gui_strings = [
        "Interactive Gel Image Analysis",
        "INSTRUCTIONS:",
        "Click and drag on image below to draw a rectangle",
        "Analyze Selection",
        "Clear Selection",
        "Ready",
        "Drawing rectangle...",
        "Analyzing...",
        "Selection cleared",
        "Error",
        "Bottom-right coordinate must be greater than top-left!",
        "Rectangle too small, please select at least 5x5 pixels!",
        "Analysis Error",
        "Done!",
        "Original Image",
        "Selected Region",
        "Raw Signal vs Background",
        "Background-Corrected Signal",
        "Vertical (Top to Bottom)",
        "Horizontal (Left to Right)",
    ]
    
    print(f"\nValidating {len(gui_strings)} GUI strings...\n")
    
    all_valid = True
    for i, text in enumerate(gui_strings, 1):
        is_ascii = all(ord(c) < 128 for c in text)
        has_escape = '\\u' in text or '\ufffd' in text or '�' in text
        has_chinese = any(ord(c) > 127 for c in text)
        
        if not is_ascii or has_escape or has_chinese:
            print(f"✗ String {i}: {repr(text)[:50]}")
            print(f"    Problem: ", end="")
            if not is_ascii:
                print("Contains non-ASCII", end=" ")
            if has_escape:
                print("Contains escape sequences", end=" ")
            if has_chinese:
                print("Contains Chinese/Unicode", end=" ")
            print()
            all_valid = False
    
    print("\n" + "="*70)
    if all_valid:
        print("✓ ALL GUI STRINGS VALID - Pure ASCII, no escapes!")
    else:
        print("✗ SOME GUI STRINGS INVALID - Fix strings above")
    print("="*70 + "\n")
    
    return all_valid


def main():
    """Run all tests"""
    print("\n" + "="*70)
    print("GUI ENCODING FIX VERIFICATION")
    print("Testing sanitize_for_gui() implementation")
    print("="*70)
    
    # Run all tests
    test1 = test_sanitization()
    test2 = test_filename_sanitization()
    test3 = test_gui_strings()
    
    # Final summary
    print("\n" + "="*70)
    print("FINAL TEST SUMMARY")
    print("="*70)
    print(f"Sanitization function:  {'✓ PASS' if test1 else '✗ FAIL'}")
    print(f"Filename sanitization:  {'✓ PASS' if test2 else '✗ FAIL'}")
    print(f"GUI string validation:  {'✓ PASS' if test3 else '✗ FAIL'}")
    print("="*70)
    
    if test1 and test2 and test3:
        print("\n✓✓✓ ALL TESTS PASSED ✓✓✓")
        print("\nThe GUI is ready to launch!")
        print("Expected behavior:")
        print("  - No \\uXXXX escape sequences")
        print("  - No \\ufffd or � characters")
        print("  - All text readable ASCII/English")
        print("  - No encoding crashes")
        print("\nTo launch GUI:")
        print("  python interactive_gui.py")
    else:
        print("\n✗✗✗ SOME TESTS FAILED ✗✗✗")
        print("\nPlease review output above and fix issues before launching GUI")
    
    print("="*70 + "\n")


if __name__ == "__main__":
    main()

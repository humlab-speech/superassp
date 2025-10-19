#!/bin/bash
# Installation script for Voice Analysis Toolbox

echo "================================================"
echo "Voice Analysis Toolbox - Python Implementation"
echo "================================================"
echo ""

# Check Python version
echo "Checking Python version..."
python3 --version

if [ $? -ne 0 ]; then
    echo "Error: Python 3 is required but not found"
    exit 1
fi

echo ""
echo "Creating virtual environment..."
python3 -m venv venv

if [ $? -ne 0 ]; then
    echo "Error: Failed to create virtual environment"
    exit 1
fi

echo "Activating virtual environment..."
source venv/bin/activate

echo ""
echo "Installing dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

if [ $? -ne 0 ]; then
    echo "Error: Failed to install dependencies"
    exit 1
fi

echo ""
echo "Installing voice_analysis package..."
pip install -e .

if [ $? -ne 0 ]; then
    echo "Error: Failed to install package"
    exit 1
fi

echo ""
echo "================================================"
echo "✓ Installation complete!"
echo "================================================"
echo ""
echo "To use the toolbox:"
echo "  1. Activate the environment: source venv/bin/activate"
echo "  2. Run analysis: python -m voice_analysis <audio_file.wav>"
echo "  3. Or use in Python:"
echo "       from voice_analysis import analyze_voice_file"
echo "       measures, F0 = analyze_voice_file('audio.wav')"
echo ""
echo "To run tests:"
echo "  python tests/test_voice_analysis.py"
echo ""
echo "To run examples:"
echo "  cd examples && python basic_usage.py"
echo ""

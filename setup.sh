#!/usr/bin/env bash
set -e

# Detect OS
OS="$(uname -s)"

# Create virtual environment
python3 -m venv venv

# Activate virtual environment
echo "Activating virtual environment..."
if [[ "$OS" == "Darwin" || "$OS" == "Linux" ]]; then
    source venv/bin/activate
elif [[ "$OS" == MINGW* || "$OS" == CYGWIN* || "$OS" == MSYS* ]]; then
    source venv/Scripts/activate
else
    echo "Please activate the virtual environment manually."
    exit 1
fi

# Upgrade pip
pip install --upgrade pip

# Install requirements
pip install -r requirements.txt

echo "Virtual environment setup complete. To activate later, use:"
if [[ "$OS" == "Darwin" || "$OS" == "Linux" ]]; then
    echo "  source venv/bin/activate"
else
    echo "  venv\\Scripts\\activate (Command Prompt) or venv\\Scripts\\Activate.ps1 (PowerShell)"
fi

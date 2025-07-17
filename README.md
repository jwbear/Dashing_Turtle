🌟 Dashing Turtle Setup & Usage Guide

Welcome! This guide will help you install and run the Dashing Turtle application, even if you have little technical background. Just follow each step carefully!
💾 Step 1: Install Docker (for the database)

Dashing Turtle needs a database to store your data. We use Docker to make this easy and automatic.
✅ Windows / macOS

    Download and install Docker Desktop: https://www.docker.com/products/docker-desktop/

    After installing, open Docker Desktop and make sure it’s running.

✅ Linux

    Install Docker Engine: Docker Engine Linux install guide

🔥 Step 2: Start the database

Open your terminal or command prompt, navigate to your project folder, and run:

docker compose up -d db

✔ This will download and start the database in the background.
✔ It will keep your data even after shutting down.
🐍 Step 3: Set up Python environment (local app)
✅ Create a virtual environment

python -m venv venv

✅ Activate it

    macOS / Linux:

source venv/bin/activate

    Windows:

venv\Scripts\activate

✅ Upgrade pip

pip install --upgrade pip

✅ Install required Python packages

pip install -r requirements.txt

✅ Install Dashing Turtle from PyPI

Dashing Turtle is also distributed via PyPI as DashingTurtle. Install it:

pip install DashingTurtle

💻 Step 4: Run the application
✅ Graphical User Interface (GUI)

dt-gui

✅ Command-Line Interface (CLI)

dt-cli

    ✅ You can choose whichever mode you prefer!

🗄️ Database notes

    Your database runs inside Docker and persists data between runs automatically.

    You can stop it any time using:

docker compose down

    To restart:

docker compose up -d db

💬 Help

Please refer to help in the command line interface for more explanation.

🔥 Data Output
Data including landscape files and figures are output to your home directory under 'DTLandscape_Output'


🔥 Sample Data for testing is available at https://github.com/jwbear/Dashing_Turtle.git



🎉 That’s it!

You’re ready to use Dashing Turtle. Enjoy exploring your data! 🚀

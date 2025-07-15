# 🌿 Dashing Turtle Quick Start Tutorial

Welcome! This guide will help you get started with sequence setup, data loading, prediction, and landscape creation.

---

## 🗺️ Overview: Visual Workflow

```plaintext
[Add Sequences]
      │
      ▼
[List Sequences & Get LIDs (Visible in Graphical Interface)]
      │
      ▼
[Load Basecall Data]
[Load Signal Data]
      │
      ▼
[Run Predictions]
      │
      ▼
[Create Landscape]

✍️ 1️⃣ Add Your Sequences

Before running predictions, you must add at least two sequences:

    One unmodified (control)

    One modified

Example command

seq -add -s GGAUCGAUCG -sec .......... -e EXP001 -n Testgoose -t 37 -t1 TypeA -t2 TypeB -r 1

Options explained:

    -s: RNA sequence (A, C, G, U)

    -sec: Secondary structure (dot-bracket notation, must match sequence length)

    -e: Experiment name (required if using secondary structure)

    -n: Sequence name

    -t: Temperature in °C (default: 37)

    -t1, -t2: Type labels (indicate modified/unmodified)

    -r: Run ID (default: 1)

🔎 2️⃣ List Sequences & Get LIDs

After adding, view your sequences and their Library IDs (LIDs):

seq -list
    52 Unmodified Example Sequence
    53 Modified Example Sequence

You will use these LIDs when loading data and running predictions.

📥 3️⃣ Load Basecall and Signal Data
Load basecall data

load -basecall -l 52 -p Sample/DMSO/Alignment

    -l: LID for this sequence

    -p: Path to basecall alignment folder

Load signal data

load -signal -l 52 -p Sample/DMSO/DMSO_fmn.txt

    -l: LID for this sequence

    -p: Path to signal file

💡 4️⃣ Run Predictions

Once both unmodified and modified data are loaded, run predictions:

predict -u 52 -l 53 -v

    -u: LID of unmodified sequence

    -l: LID of modified sequence

    -v: Include ViennaRNA base pairing probabilities

🌄 5️⃣ Create Landscape

Finally, generate a landscape visualization:

create_landscape -u 45 -l 26 -o

    -u: LID of unmodified sequence

    -l: LID of modified sequence

    -o: Optimize clusters and output dendrograms based on optimize cluster numbers

💬 Cheat Sheet
Step	Command Example
Unmodified:
Add sequence	seq -add -s AGCUAGCUA -n Test -t1 TypeA -t2 TypeA
List sequences	seq -list
Load basecall	load -basecall -l 52 -p Sample/DMSO/Alignment
Load signal	load -signal -l 52 -p /Sample/DMSO/DMSO_fmn.txt

Modified:
Add sequence	seq -add -s AGCUAGCUA -n Test -t1 TypeA -t2 TypeB
List sequences	seq -list
Load basecall	load -basecall -l 53 -p Sample/ACIM/Alignment
Load signal	load -signal -l 53 -p /Sample/DMSO/DMSO_fmn.txt

Predict	predict -u 52 -l 53 -v
Landscape	create_landscape -u 52 -l 53 -o
✅ Summary

1️⃣ Add your sequences (at least unmodified and modified)
2️⃣ List to find LIDs
3️⃣ Load data (basecall and signal)
4️⃣ Run predictions
5️⃣ Create landscape visualization

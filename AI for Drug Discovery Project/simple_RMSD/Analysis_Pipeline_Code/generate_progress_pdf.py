from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from textwrap import wrap


def draw_wrapped(c, text, x, y, max_chars=100, line_height=14):
    for line in wrap(text, width=max_chars):
        c.drawString(x, y, line)
        y -= line_height
    return y


def main():
    output_pdf = r"c:\Users\yagd9\Documents\Stony Brook AI Drug Discovery\Github\Pymol\AI for Drug Discovery Project\simple_RMSD\Analysis_Pipeline_Code\RMSD_Project_Progress_Summary_2026-07-22.pdf"

    c = canvas.Canvas(output_pdf, pagesize=letter)
    width, height = letter

    x = 0.75 * inch
    y = height - 0.9 * inch

    c.setFont("Helvetica-Bold", 16)
    c.drawString(x, y, "RMSD Workflow Progress Summary")
    y -= 20

    c.setFont("Helvetica", 11)
    c.drawString(x, y, "Date: 2026-07-22")
    y -= 22

    sections = [
        (
            "1) Scope Completed",
            [
                "Refined and simplified both RMSD tools for beginner-friendly command usage.",
                "Improved command conventions to reduce ambiguity around reference and target frame roles.",
                "Standardized documentation language for matrix orientation and mode behavior.",
            ],
        ),
        (
            "2) Script Updates",
            [
                "RMSD_per_frame_biopython_7-21-26.py: annotations expanded with clearer scenarios, fixed command ordering guidance, and matrix-vs-list mode distinctions.",
                "RMSD_calculation_function_7-21-26.py: header annotations expanded with purpose, scenarios, argument explanations, and RMSD mathematical formula.",
                "Added explicit backbone-only explanations and practical usage notes in both script docs.",
            ],
        ),
        (
            "3) README Improvements",
            [
                "README.md upgraded into a mini cookbook with detailed scenario walkthroughs.",
                "Added explicit output expectations (default filenames, list vs matrix output types).",
                "Added formula section: RMSD = sqrt((1/N) * sum(||r_i - s_i||^2)).",
                "Added troubleshooting for common Windows and dependency issues.",
                "Clarified singular vs plural reference flags:",
                "- --reference-frame => list output (Frame,RMSD_Angstrom)",
                "- --reference-frames => matrix output (references as rows, targets as columns)",
            ],
        ),
        (
            "4) Orientation and Mode Rules (Now Explicit)",
            [
                "To force Excel-style orientation (vertical references, horizontal targets), use matrix mode:",
                "--reference-frames ... --targets ...",
                "Using --reference-frame ... --targets ... does not create a matrix; it creates a simple list CSV.",
                "This was identified as a key confusion point and corrected in docs.",
            ],
        ),
        (
            "5) Commands and Execution Support Delivered",
            [
                "Provided multiple copy-pasteable PowerShell commands for 3vs1, 1vsN, NvsM, all-vs-all, and atom-filtered scenarios.",
                "Created output folder: test_outputs_for_RMSD calculations.",
                "Explained how to include self-comparisons (e.g., frame 10 vs frame 10) using matrix mode.",
            ],
        ),
        (
            "6) Environment and Error Resolution",
            [
                "Diagnosed python PATH alias issue on Windows; switched to py launcher usage.",
                "Verified Python launcher and environment (py -3.13).",
                "Resolved missing Biopython dependency and confirmed installation (Biopython 1.87).",
                "Explained harmless SyntaxWarning from backslashes in docstring examples.",
            ],
        ),
        (
            "7) Validation Status",
            [
                "Syntax/error checks run on edited Python files after modifications.",
                "No code errors reported for the updated RMSD scripts after final edits.",
            ],
        ),
        (
            "8) Current Outcome",
            [
                "The workflow now supports beginner-safe command patterns with consistent reference/target semantics.",
                "Documentation clearly separates list mode vs matrix mode and prevents axis-orientation confusion.",
                "User can now generate reproducible RMSD outputs for classroom and research reporting with minimal setup.",
            ],
        ),
    ]

    for title, bullets in sections:
        if y < 1.2 * inch:
            c.showPage()
            y = height - 0.9 * inch
            c.setFont("Helvetica", 11)

        c.setFont("Helvetica-Bold", 12)
        c.drawString(x, y, title)
        y -= 16

        c.setFont("Helvetica", 11)
        for bullet in bullets:
            prefix = "- "
            y = draw_wrapped(c, prefix + bullet, x + 8, y, max_chars=102, line_height=13)
            y -= 3

        y -= 6

    c.save()
    print(output_pdf)


if __name__ == "__main__":
    main()

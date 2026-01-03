import markdown

# Read markdown file
with open('e8-theory-of-everything/README_P_vs_NP.md', 'r', encoding='utf-8') as f:
    md_content = f.read()

# Convert to HTML
html_content = markdown.markdown(md_content, extensions=['tables', 'fenced_code'])

# Create styled HTML
styled_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>P vs NP Proof - H4 Geometric Energy Barriers</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            max-width: 1200px;
            margin: 40px auto;
            padding: 0 20px;
            line-height: 1.6;
            color: #333;
            background: #fff;
        }}
        h1 {{
            color: #1a1a1a;
            border-bottom: 3px solid #FF5722;
            padding-bottom: 10px;
            margin-top: 20px;
        }}
        h2 {{
            color: #2c3e50;
            margin-top: 40px;
            border-bottom: 2px solid #ddd;
            padding-bottom: 8px;
        }}
        h3 {{
            color: #34495e;
            margin-top: 30px;
        }}
        pre {{
            background: #f5f5f5;
            padding: 15px;
            border-left: 4px solid #FF5722;
            overflow-x: auto;
            border-radius: 4px;
            font-family: 'Courier New', Consolas, Monaco, monospace;
        }}
        code {{
            background: #f0f0f0;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', Consolas, Monaco, monospace;
            font-size: 0.9em;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
            box-shadow: 0 2px 3px rgba(0,0,0,0.1);
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background: #FF5722;
            color: white;
            font-weight: bold;
        }}
        tr:nth-child(even) {{
            background: #f9f9f9;
        }}
        tr:hover {{
            background: #f5f5f5;
        }}
        blockquote {{
            border-left: 4px solid #FF5722;
            padding-left: 20px;
            margin-left: 0;
            color: #666;
            font-style: italic;
            background: #f9f9f9;
            padding: 10px 10px 10px 20px;
            border-radius: 4px;
        }}
        a {{
            color: #FF5722;
            text-decoration: none;
        }}
        a:hover {{
            text-decoration: underline;
        }}
        hr {{
            border: none;
            border-top: 2px solid #ddd;
            margin: 40px 0;
        }}
        .header-banner {{
            background: linear-gradient(135deg, #FF5722 0%, #FFC107 100%);
            color: white;
            padding: 30px;
            border-radius: 8px;
            margin-bottom: 30px;
            text-align: center;
        }}
        .header-banner h1 {{
            border: none;
            color: white;
            margin: 0;
        }}
    </style>
</head>
<body>
    <div class="header-banner">
        <h1>🧮 P vs NP Proof</h1>
        <p>via H4 Geometric Energy Barriers</p>
    </div>
    {html_content}
</body>
</html>"""

# Write HTML file
with open('e8-theory-of-everything/README_P_vs_NP.html', 'w', encoding='utf-8') as f:
    f.write(styled_html)

print("✅ HTML file created: e8-theory-of-everything/README_P_vs_NP.html")

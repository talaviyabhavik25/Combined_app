from st_pages import Page,Section, show_pages, add_page_title
show_pages(
    [
        Page("physical_prop.py", "Physical properties", "🏠"),
        Page("comp.py", "Compressor Evaluation", "🐔"),
        Page("heater.py", "Heater Efficiency", "🐔"),
        
        # Pages after a section will be indented
        Page("line.py", "Line Sizing tool",icon="🐔"),
        # Unless you explicitly say in_section=False
        Page("pump.py", "Cent. Pump performance test",icon="🐔")
    ]
)

from st_pages import Page,Section, show_pages, add_page_title

show_pages(
    [
        Page("home_page.py", "Welcome!", "🏠"),
        
        Page("heat_exchanger.py", "Heat exchanger sizing and rating",icon="🐔"),
        Page("physical_prop.py", "Physical properties", "🐔"),
        Page("line.py", "Line Sizing tool",icon="🐔"),
        Page("comp.py", "Compressor Evaluation", "🐔"),
        Page("heater.py", "Heater Efficiency", "🐔"),
        
        # Pages after a section will be indented
        
        # Unless you explicitly say in_section=False
        Page("pump.py", "Cent. Pump performance test",icon="🐔"),
        Page('positive_pump.py', "Positive pump capacity calculations",icon="🐔")
    ]
)
from my_packages.materials_tools.materials_dataset.db_preset_scenarios import gather_entries_from_databases

ds = gather_entries_from_databases(["Mo", "Si", "B", "P"])
ds.present_as_table('test_table.txt', sort_by='Energy above hull (eV/atom)')
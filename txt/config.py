import json

config={
     "favorites"                : {"default": []},
     "unit"                     : 0,
     "save_frequency"           : 1,
     "broadening"               : "gamma_self",
     "online"                   : 1,
     "noise_filter"             : [2, 40.0],
     "manual_start_index"       : 0,
     "abundance_limit"          : 0.0,
     "legend_framealpha"        : 0.30000000000000004,
     "hold_graph_on_start"      : 1,
     "modulation_low_to_high"   : 0,
     "overlap_filter"           : [2, 50.0],
     "lower_plot"               : 5,
     "savgol_polyorder"         : 3,
     "start_index_manually"     : 0,
     "profile"                  : "Voigt",
     "full_names"               : 0,
     "uncheck_missing"          : 1,
     "vmin"                     : [1530.0, 6523.157208088715],
     "vmax"                     : [1533.0, 6535.9477124183],
     "temperature"              : [19.850000000000023, 293.0],
     "pressure"                 : [1.0, 0.9869232667160128],
     "unit_temperature"         : 1,
     "saved_tdms_column"        : "",
     "frequency"                : 100,
     "hitemp"                   : 0,
     "upper_plot"               : 5,
     "path_length"              : 11.5,
     "slot_fav"                 : "default",
     "profile_hitran_tab"       : "Voigt",
     "savgol_window_length"     : 15,
     "unit_wave"                : 0,
     "spectrum_filter"          : [0, 2.5],
     "legend_loc"               : "best",
     "unit_pressure"            : 0,
     "selections"               : { "default": []
                                    },
     "font-size"                : 4,
     "sort_molecules_by_names"  : 1,
     "slot_selection"           : "default",
     "baseline_polyorder"       : 3,
     "automatic_baseline"       : 0,
     "measure_range"            : 1.2,
     "omega_step"               : 0.01
     }


json.dump(config, open("config.txt", "w"))


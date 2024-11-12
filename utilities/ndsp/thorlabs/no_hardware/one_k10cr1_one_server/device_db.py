# replace this with your normal device db file
device_db = { # for running a test locally, i.e. with no hardware
    "core": {
        "type": "local",
        "module": "artiq.sim.devices",
        "class": "Core",
        "arguments": {}
    },
    "mains_sync": {
        "type": "local",
        "module": "artiq.sim.devices",
        "class": "Input",
        "arguments": {"name": "mains_sync"}
    }
}

# device_db["k10cr1_ndsp"] = {
#     "type": "controller",  # this tells artiq it's an ndsp
#     "sn_list": [55000741,55105674],
#     "nickname_list": ["780_QWP", "780_HWP"],
#     "host": "localhost",
#     "port": 8080
# }

device_db["k10cr1_ndsp"] = {
    "type": "controller",  # this tells artiq it's an ndsp
    "sn_list": [55000741], #,55105674],
    "nickname_list": ["780_QWP"], #, "780_HWP"],
    "host": "localhost",
    "port": 8080
}
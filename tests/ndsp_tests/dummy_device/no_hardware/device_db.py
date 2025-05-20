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

device_db["example_ndsp"] = {
    "type": "controller",  # this tells artiq it's an ndsp
    "host": "localhost",  # address
    "port": 8080
}
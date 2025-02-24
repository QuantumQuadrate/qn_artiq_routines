core_addr = "192.168.1.76"

device_db = {
    "core": {
        "type": "local",
        "module": "artiq.coredevice.core",
        "class": "Core",
        "arguments": {"host": core_addr, "ref_period": 1e-09, "target": "cortexa9"},
    },
    "core_log": {
        "type": "controller",
        "host": "::1",
        "port": 1068,
        "command": "aqctl_corelog -p {port} --bind {bind} " + core_addr
    },
    # "core_moninj": {
    #     "type": "controller",
    #     "host": "::1",
    #     "port_proxy": 1383,
    #     "port": 1384,
    #     "command": "aqctl_moninj_proxy --port-proxy {port_proxy} --port-control {port} --bind {bind} " + core_addr
    # },
    "core_cache": {
        "type": "local",
        "module": "artiq.coredevice.cache",
        "class": "CoreCache"
    },
    "core_dma": {
        "type": "local",
        "module": "artiq.coredevice.dma",
        "class": "CoreDMA"
    },

    "i2c_switch0": {
        "type": "local",
        "module": "artiq.coredevice.i2c",
        "class": "I2CSwitch",
        "arguments": {"address": 0xe0}
    },
    "i2c_switch1": {
        "type": "local",
        "module": "artiq.coredevice.i2c",
        "class": "I2CSwitch",
        "arguments": {"address": 0xe2}
    },
}

device_db["led0"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000039}
}

device_db["k10cr1_ndsp"] = {
    "type": "controller",  # this tells artiq it's an ndsp
    "sn_list": [55000741,55105674],
    "nickname_list": ["780_QWP", "780_HWP"],
    "host": "::1",  #localhost
    "port": 8080
}
#
# device_db.update({
#     "k10cr1_ndsp": {
#         "type": "controller",  # this tells artiq it's an ndsp
#         "host": "::1",
#         "port": 8080,
#         # "command": "python repository\\qn_artiq_routines\\utilities\\ndsp\\thorlabs\with_kasli\\two_k10cr1_one_server\\launcher_multi_rotor.py -p {port}"
#         "command": "python launcher_multi_rotor.py -p {port}"
#
# },
# })
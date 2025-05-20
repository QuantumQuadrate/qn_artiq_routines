import telnetlib

host = "::1"
port = 8080

tn = telnetlib.Telnet(host, port)
print("Connected to Telnet")

tn.write(b"\n")  # Send a blank line
response = tn.read_some()
print("Response:", response.decode())

tn.close()
import serial

ser= serial.Serial('COM6', 115200)

while True:
    print(ser.readline())
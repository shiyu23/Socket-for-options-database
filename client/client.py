import zmq
import json


context = zmq.Context()
socket = context.socket(zmq.SUB)
socket.connect("tcp://localhost:5555")
socket.setsockopt(zmq.SUBSCRIBE,''.encode('utf-8'))  # 接收所有消息
while True:
    response = json.loads(socket.recv().decode('utf-8'));
    print("response: %s" % response)
import zmq
import json
import re
import threading

class tcore_zmq():
    def __init__(self,APPID,SKey):
        self.context = zmq.Context()
        self.appid=APPID
        self.ServiceKey=SKey
        self.lock = threading.Lock()
        self.qsocket = None
        self.tsocket = None

#交易连线登入
    def trade_connect(self,port):
        self.lock.acquire()
        login_obj = {"Request":"LOGIN","Param":{"SystemName":self.appid, "ServiceKey":self.ServiceKey}}
        self.tsocket = self.context.socket(zmq.REQ)
        self.tsocket.connect("tcp://127.0.0.1:%s" % port)
        self.tsocket.send_string(json.dumps(login_obj))
        message = self.tsocket.recv()
        message = message[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#行情连线登入
    def quote_connect(self,port):
        self.lock.acquire()
        login_obj = {"Request":"LOGIN","Param":{"SystemName":self.appid, "ServiceKey":self.ServiceKey}}
        self.qsocket = self.context.socket(zmq.REQ)
        self.qsocket.connect("tcp://127.0.0.1:%s" % port)
        self.qsocket.send_string(json.dumps(login_obj))
        message = self.qsocket.recv()
        message = message[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#交易连线登出
    def trade_logout(self,key):
        self.lock.acquire()
        obj = {"Request":"LOGOUT","SessionKey":key}
        self.tsocket.send_string(json.dumps(obj))
        self.lock.release()
        return

#行情连线登出
    def quote_logout(self,key):
        self.lock.acquire()
        obj = {"Request":"LOGOUT","SessionKey":key}
        self.qsocket.send_string(json.dumps(obj))
        self.lock.release()
        return

#已登入资金账户
    def account_lookup(self,key):
        self.lock.acquire()
        obj = {"Request":"ACCOUNTS","SessionKey":key}
        self.tsocket.send_string(json.dumps(obj))
        message = self.tsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#查询当日委托回报
    def restore_report(self,key,page):
        self.lock.acquire()
        obj = {"Request":"RESTOREREPORT","SessionKey":key,"QryIndex":page}
        self.tsocket.send_string(json.dumps(obj))
        message = self.tsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#查询当日成交回报
    def RestoreFillReport(self,key,qryIndex):
        self.lock.acquire()
        obj = {"Request":"RESTOREFILLREPORT","SessionKey":key,"QryIndex":qryIndex}
        self.tsocket.send_string(json.dumps(obj))
        message = self.tsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#下单
    def new_order(self,key,Param):
        self.lock.acquire()
        obj = {"Request":"NEWORDER","SessionKey":key}
        obj["Param"] = Param
        self.tsocket.send_string(json.dumps(obj))
        message = self.tsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#改单
    def replace_order(self,key,Param):
        self.lock.acquire()
        obj = {"Request":"REPLACEORDER","SessionKey":key}
        obj["Param"] = Param
        self.tsocket.send_string(json.dumps(obj))
        message = self.tsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#删单
    def cancel_order(self,key,Param):
        self.lock.acquire()
        obj = {"Request":"CANCELORDER","SessionKey":key}
        obj["Param"] = Param
        self.tsocket.send_string(json.dumps(obj))
        message = self.tsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#查询资金
    def margin(self,key,AM):
        self.lock.acquire()
        obj = {"Request":"MARGINS","SessionKey":key,"AccountMask":AM}
        self.tsocket.send_string(json.dumps(obj))
        message = self.tsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#查询持仓
    def position(self,key,AM,page):
        self.lock.acquire()
        obj = {"Request":"POSITIONS","SessionKey":key,"AccountMask":AM,"QryIndex":page}
        self.tsocket.send_string(json.dumps(obj))
        message = self.tsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#订阅实时报价
    def subquote(self,key,symbol):
        self.lock.acquire()
        obj = {"Request":"SUBQUOTE","SessionKey":key}
        obj["Param"] ={"Symbol":symbol,"SubDataType":"REALTIME"}
        self.qsocket.send_string(json.dumps(obj))
        message = self.qsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#解订实时报价(每次订阅合约前，先调用解订，避免重复订阅)
    def unsubquote(self,key,symbol):
        self.lock.acquire()
        obj = {"Request":"UNSUBQUOTE","SessionKey":key}
        obj["Param"] = {"Symbol":symbol,"SubDataType":"REALTIME"}
        self.qsocket.send_string(json.dumps(obj))
        message = self.qsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#订阅实时greeks
    def subgreeks(self,key,symbol):
        self.lock.acquire()
        obj = {"Request":"SUBQUOTE","SessionKey":key}
        obj["Param"] ={"Symbol":symbol,"SubDataType":"GREEKS"}
        self.qsocket.send_string(json.dumps(obj))
        message = self.qsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#解订实时greeks(每次订阅合约前，先调用解订，避免重复订阅)
    def unsubgreeks(self,key,symbol):
        self.lock.acquire()
        obj = {"Request":"UNSUBQUOTE","SessionKey":key}
        obj["Param"] = {"Symbol":symbol,"SubDataType":"GREEKS"}
        self.qsocket.send_string(json.dumps(obj))
        message = self.qsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#订阅历史数据    
    #1：SessionKey，
    #2：合约代码，
    #3：数据周期:"TICKS","1K","DK"，
    #4: 历史数据开始时间,
    #5: 历史数据结束时间
    def sub_history(self,key,Param):
        self.lock.acquire()
        obj = {"Request":"SUBQUOTE","SessionKey":key}
        obj["Param"] = Param
        self.qsocket.send_string(json.dumps(obj))
        message = self.qsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data 

#解订历史数据（遗弃，不再使用）
    #1：SessionKey，
    #2：合约代码，
    #3：数据周期"TICKS","1K","DK"，
    #4: 历史数据开始时间,
    #5: 历史数据结束时间   
    def un_subhistory(self,key,symbol,type,starttime,endtime):
        self.lock.acquire()
        obj = {"Request":"UNSUBQUOTE","SessionKey":key}
        obj["Param"] ={"Symbol": symbol,"SubDataType":type,"StartTime" :starttime,"EndTime" :endtime}
        self.qsocket.send_string(json.dumps(obj))
        message = self.qsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#分页获取订阅的历史数据
    def get_history(self,key,Param):
        self.lock.acquire()
        obj = {"Request":"GETHISDATA","SessionKey":key}
        obj["Param"] = Param
        self.qsocket.send_string(json.dumps(obj))
        message = (self.qsocket.recv()[:-1]).decode("utf-8")
        index =  re.search(":",message).span()[1]  # filter 
        message = message[index:]
        message = json.loads(message)
        self.lock.release()
        return message

#查询合约信息
    def QueryInstrumentInfo(self, key, sym):
        self.lock.acquire()
        obj = {"Request" : "QUERYINSTRUMENTINFO" , "SessionKey" : key , "Symbol" : sym}
        self.qsocket.send_string(json.dumps(obj))
        message = self.qsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#查询对应类型的所有合约
    #"Type":
    #期货：Future
    #期权：Options
    #证券：Stock
    def QueryAllInstrumentInfo(self, key, type):
        self.lock.acquire()
        obj = {"Request": "QUERYALLINSTRUMENT", "SessionKey": key, "Type": type}
        self.qsocket.send_string(json.dumps(obj))
        message = self.qsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#交易连线心跳（在收到"PING"消息时调用）
    def TradePong(self,key, id = ""):
        if self.tsocket == None:
            return

        self.lock.acquire()
        obj = {"Request":"PONG","SessionKey":key, "ID":id}
        self.tsocket.send_string(json.dumps(obj))
        message = self.tsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

#行情连线心跳（在收到"PING"消息时调用）
    def QuotePong(self,key, id = ""):
        if self.qsocket == None:
            return

        self.lock.acquire()
        obj = {"Request":"PONG","SessionKey":key, "ID":id}
        self.qsocket.send_string(json.dumps(obj))
        message = self.qsocket.recv()[:-1]
        data = json.loads(message)
        self.lock.release()
        return data

class KeepAliveHelper():
    def __init__(self, sub_port, session, objZMQ):
        threading.Thread(target = self.ThreadProcess, args=(sub_port, session, objZMQ)).start()
        self.IsTerminal = False

    def Close(self):
        self.IsTerminal = True

    def ThreadProcess(self, sub_port, session, objZMQ):
        socket_sub = zmq.Context().socket(zmq.SUB)
        socket_sub.connect("tcp://127.0.0.1:%s" % sub_port)
        socket_sub.setsockopt_string(zmq.SUBSCRIBE,"")
        while True:
            message = (socket_sub.recv()[:-1]).decode("utf-8")
            findText = re.search("{\"DataType\":\"PING\"}",message)

            if findText == None:
                continue

            if self.IsTerminal:
                return

            objZMQ.QuotePong(session, "TC")
            objZMQ.TradePong(session, "TC")
import zmq
import time
import sys

sys.coinit_flags = 0
from tcoreapi_mq import *
import re

import math
import numpy as np
from scipy.interpolate import CubicSpline
import calendar
from enum import Enum
import os


class Maturity(Enum):
    M1 = 1; M2 = 2; M3 = 3; Q1 = 4; Q2 = 5; Q3 = 6

class StockType(Enum):
    etf50 = 1; h300 = 2; gz300 = 3; s300 = 4

class FutureType(Enum):
    IF = 1; IH = 2

class OptionType(Enum):
    C = 1; P = 2

holiday = (calendar.datetime.date(2021, 1, 1),
                        calendar.datetime.date(2021, 2, 11),
                        calendar.datetime.date(2021, 2, 12),
                        calendar.datetime.date(2021, 2, 15),
                        calendar.datetime.date(2021, 2, 16),
                        calendar.datetime.date(2021, 2, 17),
                        calendar.datetime.date(2021, 4, 5),
                        calendar.datetime.date(2021, 5, 3),
                        calendar.datetime.date(2021, 5, 4),
                        calendar.datetime.date(2021, 5, 5),
                        calendar.datetime.date(2021, 6, 14),
                        calendar.datetime.date(2021, 9, 20),
                        calendar.datetime.date(2021, 9, 21),
                        calendar.datetime.date(2021, 10, 1),
                        calendar.datetime.date(2021, 10, 4),
                        calendar.datetime.date(2021, 10, 5),
                        calendar.datetime.date(2021, 10, 6),
                        calendar.datetime.date(2021, 10, 7),

                        calendar.datetime.date(2022, 1, 3),
                        calendar.datetime.date(2022, 1, 31),
                        calendar.datetime.date(2022, 2, 1),
                        calendar.datetime.date(2022, 2, 2),
                        calendar.datetime.date(2022, 2, 3),
                        calendar.datetime.date(2022, 2, 4),
                        calendar.datetime.date(2022, 4, 4),
                        calendar.datetime.date(2022, 4, 5),
                        calendar.datetime.date(2022, 5, 2),
                        calendar.datetime.date(2022, 5, 3),
                        calendar.datetime.date(2022, 5, 4),
                        calendar.datetime.date(2022, 5, 5),
                        calendar.datetime.date(2022, 6, 3),
                        calendar.datetime.date(2022, 9, 12),
                        calendar.datetime.date(2022, 10, 3),
                        calendar.datetime.date(2022, 10, 4),
                        calendar.datetime.date(2022, 10, 5),
                        calendar.datetime.date(2022, 10, 6),
                        calendar.datetime.date(2022, 10, 7),
    ) # 2021 + 2022


def cdf(x: float):

    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911
    sign = 1
    if x < 0:
        sign = -1
    x = math.fabs(x) / math.sqrt(2.0);
    t = 1.0 / (1.0 + p * x);
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(- x * x)
    return 0.5 * (1.0 + sign * y)


def pdf(x: float):

    pi = 3.141592653589793
    return 1 / (math.sqrt(2 * pi)) * math.exp(- x * x / 2)


def BS(oty: OptionType, K: float, T: float, S: float, sigma: float):

    sigmaSqrtT = sigma * math.sqrt(T)
    d1 = math.log(S / K) / sigmaSqrtT + 0.5 * sigmaSqrtT
    d2 = d1 - sigmaSqrtT
    if oty == OptionType.C:
        return S * cdf(d1) - K * cdf(d2)
    else:
        return K * cdf(-d2) - S * cdf(-d1)


class OptionInfo:

    def __init__(self, sty: StockType, mat: Maturity, oty: OptionType, K: float):
        self.sty = sty
        self.mat = mat
        self.oty = oty
        self.K = K
        self.P: float
        self.P_yc: float
        self.ask: float = ''
        self.bid: float = ''
        self.T: float
        self.S: float
        self.cb: bool = False
        self._iv: float
        self.settlement_price: float
        self.deposit: float
        self.yc_master_contract: str = ''

    def midbidaskspread(self):
        if '' not in [self.ask, self.bid]:
            return (self.ask + self.bid)/2
        else:
            return ''

    def iv(self):
        a = 0.0001; b = 3; NTRY = 20; FACTOR = 1.6; S = self.S; T = self.T; K = self.K; P = self.midbidaskspread(); oty = self.oty
        f1 = BS(oty, K, T, S, a) - P; f2 = BS(oty, K, T, S, b) - P
        # rfbisect
        tol = 1e-6
        while (b - a) > tol and NTRY > 0:
            NTRY -= 1
            c = (a + b) / 2.0
            if abs(BS(oty, K, T, S, c) - P) < tol:
                return c
            else:
                if (BS(oty, K, T, S, a) - P) * (BS(oty, K, T, S, c) - P) < 0:
                    b = c
                else:
                    a = c
        return c

    def delta(self):
        iv = self._iv; S = self.S; T = self.T
        if self.oty == OptionType.C:
            return cdf(math.log(S / self.K) / (iv * math.sqrt(T)) + 0.5 * iv * math.sqrt(T))
        else:
            return cdf(math.log(S / self.K) / (iv * math.sqrt(T)) + 0.5 * iv * math.sqrt(T)) - 1

    def gamma(self):
        iv = self._iv; S = self.S; T = self.T
        return pdf(math.log(S / self.K) / (iv * math.sqrt(T)) + 0.5 * iv * math.sqrt(T)) / S / iv / math.sqrt(T)

    def vega(self):
        iv = self._iv; S = self.S; T = self.T
        return S * math.sqrt(T) * pdf(math.log(S / self.K) / (iv * math.sqrt(T)) + 0.5 * iv * math.sqrt(T))

    def theta(self):
        iv = self._iv; S = self.S; T = self.T
        return - S * pdf(math.log(S / self.K) / (iv * math.sqrt(T)) + 0.5 * iv * math.sqrt(T)) * iv / 2 / math.sqrt(T)

    def _deposit(self, ul_yc: float):
        if ul_yc == '':
            return ''

        if self.sty == StockType.gz300:
            if self.oty == OptionType.C:
                self.deposit = self.settlement_price * 100 + max(ul_yc * 100 * 0.1 - max(self.K - ul_yc, 0) * 100, 0.5 * ul_yc * 100 * 0.1)
            elif self.oty == OptionType.P:
                self.deposit = self.settlement_price * 100 + max(ul_yc * 100 * 0.1 - max(ul_yc - self.K, 0) * 100, 0.5 * self.K * 100 * 0.1)
        elif self.sty in [StockType.etf50, StockType.h300, StockType.s300]:
            if self.oty == OptionType.C:
                self.deposit = (self.settlement_price + max(0.12 * ul_yc - max(self.K - ul_yc, 0), 0.07 * ul_yc)) * 10000
            elif self.oty == OptionType.P:
                self.deposit = min(self.settlement_price + max(0.12 * ul_yc - max(ul_yc - self.K, 0), 0.07 * self.K), self.K) * 10000


class OptData:

    def __init__(self, sty: StockType):
        self.sty = sty
        self.Mat_to_2005 = {}
        self._2005_to_Mat = {}
        self.T = {}
        self.initT = {}
        self.S = {}
        self.k0 = {}
        self.posi = {}
        self.OptionList = {}
        self.ul_yc: float = ''
        if sty == StockType.gz300:
            self.cm = 100
            self.mc = 0.2
            self.matlist = [Maturity.M1, Maturity.M2, Maturity.M3, Maturity.Q1, Maturity.Q2, Maturity.Q3]
        elif sty in [StockType.etf50, StockType.h300, StockType.s300]:
            self.cm = 10000
            self.mc = 0.0001
            self.matlist = [Maturity.M1, Maturity.M2, Maturity.Q1, Maturity.Q2]
        for mat in self.matlist:
            self.S[mat] = ''; self.k0[mat] = ''; self.posi[mat] = ''
        self.k_list = {}
        self.getMat()
        
    def getMat(self):
        for mat in self.matlist:
            self.Mat_to_2005[mat] = Mat['contract_format'][self.sty][mat]
            self._2005_to_Mat[self.Mat_to_2005[mat]] = mat

        def num_weekend(date1: calendar.datetime.date, date2: calendar.datetime.date):
            num = 0
            oneday = calendar.datetime.timedelta(days = 1)
            date = calendar.datetime.date(date1.year, date1.month, date1.day)
            while date != date2:
                if date.weekday() == 5 or date.weekday() == 6 or date in holiday:
                    num += 1
                date += oneday
            return num

        c = calendar.Calendar(firstweekday=calendar.SUNDAY)
        localtime = time.localtime()
        year = localtime.tm_year; month = localtime.tm_mon; mday = localtime.tm_mday
        currentDate = calendar.datetime.date(year, month, mday)

        for mat in self.matlist:
            self.T[mat] = ((Mat['calendar'][self.sty][mat] - currentDate).days - num_weekend(currentDate, Mat['calendar'][self.sty][mat]))/244
            self.initT[mat] = self.T[mat]

    def subscribe_init(self, mat: Maturity):
        QuoteID_addin = []
        for id in QuoteID:
            if self.sty == StockType.gz300 and id[11:13] == 'IO' and id[16:20] == self.Mat_to_2005[mat]:
                QuoteID_addin.append(id)
            elif self.sty == StockType.etf50 and id[9:15] == '510050' and self.Mat_to_2005[mat] in [id[18:22], id[19:23]]:
                QuoteID_addin.append(id)
            elif self.sty == StockType.h300 and id[9:15] == '510300' and self.Mat_to_2005[mat] in [id[18:22], id[19:23]]:
                QuoteID_addin.append(id)
            elif self.sty == StockType.s300 and id[10:16] == '159919' and self.Mat_to_2005[mat] in [id[19:23], id[20:24]]:
                QuoteID_addin.append(id)

        self.OptionList[mat] = []
        QuoteID_addin_C_K = []
        for id in QuoteID_addin:
            if '.C.' in id:
                QuoteID_addin_C_K.append(float(id[last_C_P(id):]))
        QuoteID_addin_C_K.sort()
        for k in QuoteID_addin_C_K:
            self.OptionList[mat].append([OptionInfo(self.sty, mat, OptionType.C, k), OptionInfo(self.sty, mat, OptionType.P, k)])

        self.k_list[mat] = QuoteID_addin_C_K

    def S_k0_posi(self, mat: Maturity):
        optlist = self.OptionList[mat]
        n = len(optlist)
        future = [optlist[i][0].midbidaskspread() - optlist[i][1].midbidaskspread() + optlist[i][0].K for i in range(n) if '' not in [optlist[i][0].midbidaskspread(), optlist[i][1].midbidaskspread()] and 'A' not in optlist[i][0].yc_master_contract]
        future.sort()
        if future[1:-1] == []:
            return
        avg = np.mean(future[1:-1])
        self.S[mat] = avg
        self.posi[mat] = np.argmin(abs(np.array(self.k_list[mat]) - avg))
        self.k0[mat] = optlist[self.posi[mat]][0].K

    def vix(self, mat: Maturity):
        k_list_copy = self.k_list[mat].copy()
        k1 = k_list_copy[np.argmin(abs(np.array(k_list_copy) - self.S[mat]))]
        k_list_copy.remove(k1)
        k2 = k_list_copy[np.argmin(abs(np.array(k_list_copy) - self.S[mat]))]
        k_list_copy.remove(k2)
        k3 = k_list_copy[np.argmin(abs(np.array(k_list_copy) - self.S[mat]))]
        k_list_copy.remove(k3)

        [k1, k2, k3] = sorted([k1, k2, k3])

        opt1 = self.OptionList[mat][self.k_list[mat].index(k1)]
        opt2 = self.OptionList[mat][self.k_list[mat].index(k2)]
        opt3 = self.OptionList[mat][self.k_list[mat].index(k3)]
        x = [opt1[0].K, opt2[0].K, opt3[0].K]
        y = [(opt1[0].iv() + opt1[1].iv()) / 2, (opt2[0].iv() + opt2[1].iv()) / 2, (opt3[0].iv() + opt3[1].iv()) / 2]
        cs = CubicSpline(x, y)
        return cs(self.S[mat])

    def skew_same_T(self, mat: Maturity):
        optlist = self.OptionList[mat]
        n = len(optlist)
        f0 = self.S[mat]
        k0 = self.k0[mat]
        p1 = - (1 + math.log(f0 / k0) - f0 / k0); p2 = 2 * math.log(k0 / f0) * (f0 / k0 - 1) + 1/2 * math.log(k0 / f0) ** 2; p3 = 3 * math.log(k0 / f0) ** 2 * (1/3 * math.log(k0 / f0) - 1 + f0 / k0)
        for i in range(n):
            if optlist[i][0].K <= f0:
                if i == 0:
                    p1 += - 1 / (optlist[i][1].K) ** 2 * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i][1].K)
                    p2 += 2 / (optlist[i][1].K) ** 2 * (1 - math.log(optlist[i][1].K / f0)) * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i][1].K)
                    p3 += 3 / (optlist[i][1].K) ** 2 * (2 * math.log(optlist[i][1].K / f0) - math.log(optlist[i][1].K / f0) ** 2) * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i][1].K)
                elif i == (n-1):
                    p1 += - 1 / (optlist[i][1].K) ** 2 * optlist[i][1].midbidaskspread() * (optlist[i][1].K - optlist[i - 1][1].K)
                    p2 += 2 / (optlist[i][1].K) ** 2 * (1 - math.log(optlist[i][1].K / f0)) * optlist[i][1].midbidaskspread() * (optlist[i][1].K - optlist[i - 1][1].K)
                    p3 += 3 / (optlist[i][1].K) ** 2 * (2 * math.log(optlist[i][1].K / f0) - math.log(optlist[i][1].K / f0) ** 2) * optlist[i][1].midbidaskspread() * (optlist[i][1].K - optlist[i - 1][1].K)
                else:
                    p1 += - 1 / (optlist[i][1].K) ** 2 * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i - 1][1].K) / 2
                    p2 += 2 / (optlist[i][1].K) ** 2 * (1 - math.log(optlist[i][1].K / f0)) * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i - 1][1].K) / 2
                    p3 += 3 / (optlist[i][1].K) ** 2 * (2 * math.log(optlist[i][1].K / f0) - math.log(optlist[i][1].K / f0) ** 2) * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i - 1][1].K) / 2
            elif optlist[i][0].K >= f0:
                if i == (n-1):
                    p1 += - 1 / (optlist[i][0].K) ** 2 * optlist[i][0].midbidaskspread() * (optlist[i][0].K - optlist[i - 1][0].K)
                    p2 += 2 / (optlist[i][0].K) ** 2 * (1 - math.log(optlist[i][0].K / f0)) * optlist[i][0].midbidaskspread() * (optlist[i][0].K - optlist[i - 1][0].K)
                    p3 += 3 / (optlist[i][0].K) ** 2 * (2 * math.log(optlist[i][0].K / f0) - math.log(optlist[i][0].K / f0) ** 2) * optlist[i][0].midbidaskspread() * (optlist[i][0].K - optlist[i - 1][0].K)
                elif i == 0:
                    p1 += - 1 / (optlist[i][0].K) ** 2 * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i][0].K)
                    p2 += 2 / (optlist[i][0].K) ** 2 * (1 - math.log(optlist[i][0].K / f0)) * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i][0].K)
                    p3 += 3 / (optlist[i][0].K) ** 2 * (2 * math.log(optlist[i][0].K / f0) - math.log(optlist[i][0].K / f0) ** 2) * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i][0].K)
                else:
                    p1 += - 1 / (optlist[i][0].K) ** 2 * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i - 1][0].K) / 2
                    p2 += 2 / (optlist[i][0].K) ** 2 * (1 - math.log(optlist[i][0].K / f0)) * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i - 1][0].K) / 2
                    p3 += 3 / (optlist[i][0].K) ** 2 * (2 * math.log(optlist[i][0].K / f0) - math.log(optlist[i][0].K / f0) ** 2) * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i - 1][0].K) / 2
        if p2 - p1 ** 2 > 0:
            return (p3 - 3 * p1 * p2 + 2 * p1 ** 3) / math.sqrt((p2 - p1 ** 2) ** 3)
        else:
            return 0

    def skew_same_T_partial(self, mat: Maturity, partial: int):
        optlist = self.OptionList[mat]
        n = len(optlist)
        f0 = self.S[mat]
        k0 = self.k0[mat]
        cen = self.posi[mat]
        p1 = - (1 + math.log(f0 / k0) - f0 / k0); p2 = 2 * math.log(k0 / f0) * (f0 / k0 - 1) + 1/2 * math.log(k0 / f0) ** 2; p3 = 3 * math.log(k0 / f0) ** 2 * (1/3 * math.log(k0 / f0) - 1 + f0 / k0)
        for i in range(max(cen - partial, 0), min(cen + partial + 1, n)):
            if optlist[i][0].K <= f0:
                if i == 0:
                    p1 += - 1 / (optlist[i][1].K) ** 2 * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i][1].K)
                    p2 += 2 / (optlist[i][1].K) ** 2 * (1 - math.log(optlist[i][1].K / f0)) * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i][1].K)
                    p3 += 3 / (optlist[i][1].K) ** 2 * (2 * math.log(optlist[i][1].K / f0) - math.log(optlist[i][1].K / f0) ** 2) * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i][1].K)
                elif i == (n-1):
                    p1 += - 1 / (optlist[i][1].K) ** 2 * optlist[i][1].midbidaskspread() * (optlist[i][1].K - optlist[i - 1][1].K)
                    p2 += 2 / (optlist[i][1].K) ** 2 * (1 - math.log(optlist[i][1].K / f0)) * optlist[i][1].midbidaskspread() * (optlist[i][1].K - optlist[i - 1][1].K)
                    p3 += 3 / (optlist[i][1].K) ** 2 * (2 * math.log(optlist[i][1].K / f0) - math.log(optlist[i][1].K / f0) ** 2) * optlist[i][1].midbidaskspread() * (optlist[i][1].K - optlist[i - 1][1].K)
                else:
                    p1 += - 1 / (optlist[i][1].K) ** 2 * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i - 1][1].K) / 2
                    p2 += 2 / (optlist[i][1].K) ** 2 * (1 - math.log(optlist[i][1].K / f0)) * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i - 1][1].K) / 2
                    p3 += 3 / (optlist[i][1].K) ** 2 * (2 * math.log(optlist[i][1].K / f0) - math.log(optlist[i][1].K / f0) ** 2) * optlist[i][1].midbidaskspread() * (optlist[i + 1][1].K - optlist[i - 1][1].K) / 2
            elif optlist[i][0].K >= f0:
                if i == (n-1):
                    p1 += - 1 / (optlist[i][0].K) ** 2 * optlist[i][0].midbidaskspread() * (optlist[i][0].K - optlist[i - 1][0].K)
                    p2 += 2 / (optlist[i][0].K) ** 2 * (1 - math.log(optlist[i][0].K / f0)) * optlist[i][0].midbidaskspread() * (optlist[i][0].K - optlist[i - 1][0].K)
                    p3 += 3 / (optlist[i][0].K) ** 2 * (2 * math.log(optlist[i][0].K / f0) - math.log(optlist[i][0].K / f0) ** 2) * optlist[i][0].midbidaskspread() * (optlist[i][0].K - optlist[i - 1][0].K)
                elif i == 0:
                    p1 += - 1 / (optlist[i][0].K) ** 2 * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i][0].K)
                    p2 += 2 / (optlist[i][0].K) ** 2 * (1 - math.log(optlist[i][0].K / f0)) * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i][0].K)
                    p3 += 3 / (optlist[i][0].K) ** 2 * (2 * math.log(optlist[i][0].K / f0) - math.log(optlist[i][0].K / f0) ** 2) * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i][0].K)
                else:
                    p1 += - 1 / (optlist[i][0].K) ** 2 * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i - 1][0].K) / 2
                    p2 += 2 / (optlist[i][0].K) ** 2 * (1 - math.log(optlist[i][0].K / f0)) * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i - 1][0].K) / 2
                    p3 += 3 / (optlist[i][0].K) ** 2 * (2 * math.log(optlist[i][0].K / f0) - math.log(optlist[i][0].K / f0) ** 2) * optlist[i][0].midbidaskspread() * (optlist[i + 1][0].K - optlist[i - 1][0].K) / 2
        if p2 - p1 ** 2 > 0:
            return (p3 - 3 * p1 * p2 + 2 * p1 ** 3) / math.sqrt((p2 - p1 ** 2) ** 3)
        else:
            return 0

    def skew(self, mat1: Maturity, mat2: Maturity):
        w = (self.T[mat2] - 1/12) / (self.T[mat2] - self.T[mat1])
        s1 = self.skew_same_T(mat1)
        s2 = self.skew_same_T(mat2)
        return 100 - 10 * (w * s1 + (1 - w) * s2)

    def set_deposit(self, opt: OptionInfo):
        opt._deposit(self.ul_yc)

class FutureData:

    def __init__(self, fty: FutureType):
        self.fty = fty
        self.Mat_to_2005 = {}
        self._2005_to_Mat = {}
        self.initT = {}
        self.cm = 300
        self.matlist = [Maturity.M1, Maturity.M2, Maturity.Q1, Maturity.Q2]
        self.getMat()
        
    def getMat(self):
        for mat in self.matlist:
            self.Mat_to_2005[mat] = Mat['contract_format'][self.fty][mat]
            self._2005_to_Mat[self.Mat_to_2005[mat]] = mat

        def num_weekend(date1: calendar.datetime.date, date2: calendar.datetime.date):
            num = 0
            oneday = calendar.datetime.timedelta(days = 1)
            date = calendar.datetime.date(date1.year, date1.month, date1.day)
            while date != date2:
                if date.weekday() == 5 or date.weekday() == 6 or date in holiday:
                    num += 1
                date += oneday
            return num

        c = calendar.Calendar(firstweekday=calendar.SUNDAY)
        localtime = time.localtime()
        year = localtime.tm_year; month = localtime.tm_mon; mday = localtime.tm_mday
        currentDate = calendar.datetime.date(year, month, mday)

        for mat in self.matlist:
            self.initT[mat] = ((Mat['calendar'][self.fty][mat] - currentDate).days - num_weekend(currentDate, Mat['calendar'][self.fty][mat]))/244


class Server:

    def __init__(self, parent=None):
        context = zmq.Context()
        self.socket = context.socket(zmq.PUB)
        self.socket.bind("tcp://*:5555")
        self.tydict = {StockType.etf50: 'ETF50', StockType.gz300: 'GZ300', StockType.h300: 'H300', StockType.s300: 'S300', FutureType.IF: 'IF', FutureType.IH: 'IH'}
        self.matdict = {Maturity.M1: 'M1', Maturity.M2: 'M2', Maturity.M3: 'M3', Maturity.Q1: 'Q1', Maturity.Q2: 'Q2', Maturity.Q3: 'Q3'}
        self.otydict = {OptionType.C: 'C', OptionType.P: 'P'}

    def send(self, msg):
        if msg == 'end':
            os._exit(0)
        self.socket.send(bytes(json.dumps(msg).encode('utf-8')))

    def get_and_cal(self, quote):

        QuoteID = quote["Symbol"]
        TradingPrice = quote['TradingPrice'] if quote['TradingPrice'] == '' else float(quote['TradingPrice'])
        Bid = quote['Bid'] if quote['Bid'] == '' else float(quote['Bid'])
        Ask = quote['Ask'] if quote['Ask'] == '' else float(quote['Ask'])
        YClosedPrice = quote['YClosedPrice'] if quote['YClosedPrice'] == '' else float(quote['YClosedPrice'])
        HighPrice = quote['HighPrice'] if quote['HighPrice'] == '' else float(quote['HighPrice'])
        LowPrice = quote['LowPrice'] if quote['LowPrice'] == '' else float(quote['LowPrice'])

        t = time.localtime()

        # option
        if QuoteID[3] == 'O':

            opt = name_to_data(QuoteID)

            sty = opt.sty
            mat = opt.mat

            # update OptionList
            opt.P = TradingPrice
            opt.bid = Bid
            opt.ask = Ask
            opt.yc_master_contract = QuoteID
            opt.settlement_price = quote['ReferencePrice'] if quote['ReferencePrice'] == '' else float(quote['ReferencePrice'])
            # update S, k0, posi
            data_opt[sty].S_k0_posi(mat)
            opt.S = data_opt[sty].S[mat]
            # update time
            data_opt[sty].T[mat] = data_opt[sty].initT[mat] + ((15 - t.tm_hour - 1 - 1.5 * (t.tm_hour < 12)) * 60 * 60 + (60 - t.tm_min -1) * 60 + (60 - t.tm_sec) + 120) / (60 * 60 * 4 + 120) / 244
            opt.T = data_opt[sty].T[mat]
            # cb
            try:
                if opt.cb == False and float(quote["Bid"]) == float(quote["Ask"]):
                    opt.cb = True
                if opt.cb == True and not float(quote["Bid"]) == float(quote["Ask"]):
                    opt.cb = False
            except:
                pass

            S = opt.S
            MID = (Ask + Bid) / 2 if '' not in [Ask, Bid] else ''
            IV = ''
            if '' not in [S, MID]:
                opt._iv = opt.iv()
                IV = opt._iv

            msg = {'Type': 'Option',
                            'CONTRACT': QuoteID,
                            'TIME': time.strftime('%Y/%m/%d %H:%M:%S', t),
                            'T': opt.T,
                            'STY': self.tydict[sty],
                            'MAT': self.matdict[mat],
                            'MAT_DATE': Mat['calendar'][sty][mat].strftime('%Y/%m/%d'),
                            'OTY': self.otydict[opt.oty],
                            'K': opt.K,
                            'P': TradingPrice,
                            'P_YC': YClosedPrice,
                            'ASK': Ask,
                            'BID': Bid,
                            'MID': MID,
                            'P_HIGHEST': HighPrice,
                            'P_LOWEST': LowPrice,
                            'S': S,
                            'SETTLEMENT_PRICE': opt.settlement_price,
                            'DEPOSIT': opt._deposit(data_opt[sty].ul_yc),
                            'IV': IV,
                            'DELTA': opt.delta() if not IV == '' else '',
                            'GAMMA': opt.gamma() if not IV == '' else '',
                            'VEGA': opt.vega() if not IV == '' else '',
                            'THETA': opt.theta() if not IV == '' else '',
                   }
            self.send(msg)


        # future
        elif QuoteID[3] == 'F':

            if 'IF' in QuoteID:
                fty = FutureType.IF
            elif 'IH' in QuoteID:
                fty = FutureType.IH
            else:
                return

            mat = data_opt[fty]._2005_to_Mat[QuoteID[-4:]]
            SETTLEMENT_PRICE = quote['ReferencePrice'] if quote['ReferencePrice'] == '' else float(quote['ReferencePrice'])

            msg = {'Type': 'Future',
                            'CONTRACT': QuoteID,
                            'TIME': time.strftime('%Y/%m/%d %H:%M:%S', t),
                            'T': data_opt[fty].initT[mat] + ((15 - t.tm_hour - 1 - 1.5 * (t.tm_hour < 12)) * 60 * 60 + (60 - t.tm_min -1) * 60 + (60 - t.tm_sec) + 120) / (60 * 60 * 4 + 120) / 244,
                            'FTY': self.tydict[fty],
                            'MAT': self.matdict[mat],
                            'MAT_DATE': Mat['calendar'][fty][mat].strftime('%Y/%m/%d'),
                            'P': TradingPrice,
                            'P_YC': YClosedPrice,
                            'ASK': Ask,
                            'BID': Bid,
                            'MID': (Ask + Bid) / 2 if '' not in [Ask, Bid] else '',
                            'P_HIGHEST': HighPrice,
                            'P_LOWEST': LowPrice,
                            'SETTLEMENT_PRICE': SETTLEMENT_PRICE,
                            'DEPOSIT': SETTLEMENT_PRICE * 300 * 0.12 if not SETTLEMENT_PRICE == '' else '',
                   }
            self.send(msg)


        # underlying
        elif QuoteID[3] == 'S':

            if '510050' in QuoteID:
                sty = StockType.etf50
            elif '510300' in QuoteID:
                sty = StockType.h300
            elif '159919' in QuoteID:
                sty = StockType.s300
            elif '000300' in QuoteID:
                sty = StockType.gz300
            else:
                return

            data_opt[sty].ul_yc = YClosedPrice

            msg = {'Type': 'Stock',
                            'CONTRACT': QuoteID,
                            'TIME': time.strftime('%Y/%m/%d %H:%M:%S', t),
                            'STY': self.tydict[sty],
                            'P': TradingPrice,
                            'P_YC': YClosedPrice,
                            'ASK': Ask,
                            'BID': Bid,
                            'MID': (Ask + Bid) / 2 if '' not in [Ask, Bid] else '',
                            'P_HIGHEST': HighPrice,
                            'P_LOWEST': LowPrice,
                   }
            self.send(msg)


def last_C_P(string: str):
    num = len(string) - 1
    while (string[num] != 'C' and string[num] != 'P'):
        num -= 1
    return num + 2


def name_to_data(yc_master_contract: str):

    if 'TC.O.SSE.510050' in yc_master_contract:
        sty = StockType.etf50
        if yc_master_contract[15] == 'A':
            mat = data_opt[sty]._2005_to_Mat[yc_master_contract[19 : 23]]
        else:
            mat = data_opt[sty]._2005_to_Mat[yc_master_contract[18 : 22]]
    elif 'TC.O.SSE.510300' in yc_master_contract:
        sty = StockType.h300
        if yc_master_contract[15] == 'A':
            mat = data_opt[sty]._2005_to_Mat[yc_master_contract[19 : 23]]
        else:
            mat = data_opt[sty]._2005_to_Mat[yc_master_contract[18 : 22]]
    elif 'TC.O.SZSE.159919' in yc_master_contract:
        sty = StockType.s300
        if yc_master_contract[16] == 'A':
            mat = data_opt[sty]._2005_to_Mat[yc_master_contract[20 : 24]]
        else:
            mat = data_opt[sty]._2005_to_Mat[yc_master_contract[19 : 23]]
    elif 'TC.O.CFFEX.IO' in yc_master_contract:
        sty = StockType.gz300
        mat = data_opt[sty]._2005_to_Mat[yc_master_contract[16 : 20]]

    position = data_opt[sty].k_list[mat].index(float(yc_master_contract[last_C_P(yc_master_contract) : ]))
    if '.C.' in yc_master_contract:
        se = 0
    elif '.P.' in yc_master_contract:
        se = 1

    return data_opt[sty].OptionList[mat][position][se]


def sub_all_options():

    global g_QuoteZMQ
    global g_QuoteSession
    global q_data
    g_QuoteZMQ = None
    g_QuoteSession = ''
    q_data = None

    global QuoteID
    QuoteID = []

    global Mat
    Mat = {'calendar': {}, 'contract_format': {}}
    for format in ['calendar', 'contract_format']:
        for ty in [StockType.gz300, StockType.etf50, StockType.h300, StockType.s300, FutureType.IF, FutureType.IH]:
            Mat[format][ty] = []

    global data_opt
    data_opt = {}


    g_QuoteZMQ = tcore_zmq("ZMQ","8076c9867a372d2a9a814ae710c256e2")
    q_data = g_QuoteZMQ.quote_connect("51878") # 方正 51909、51522，公版 51878

    if q_data["Success"] != "OK":
        print("[quote]connection failed")
        return
    g_QuoteSession = q_data["SessionKey"]

    data = g_QuoteZMQ.QueryAllInstrumentInfo(g_QuoteSession, "Options")
    for i in range(len(data['Instruments']["Node"])):
        if data['Instruments']["Node"][i]['ENG'] == 'SSE(O)':
            for mat_classification in data['Instruments']["Node"][i]["Node"][0]["Node"][-4 : ]:
                for z in range(2):
                    QuoteID += mat_classification["Node"][z]['Contracts'] # etf50; z =1 for call; z=2 for put
                    Mat['calendar'][StockType.etf50] += [calendar.datetime.date(int(x[0:4]), int(x[4:6]), int(x[-2:])) for x in mat_classification["Node"][z]['ExpirationDate']]
                    Mat['contract_format'][StockType.etf50] += [x[2:6] for x in mat_classification["Node"][z]['ExpirationDate']]
            for mat_classification in data['Instruments']["Node"][i]["Node"][1]["Node"][-4 : ]:
                for z in range(2):
                    QuoteID += mat_classification["Node"][z]['Contracts'] # h300
                    Mat['calendar'][StockType.h300] += [calendar.datetime.date(int(x[0:4]), int(x[4:6]), int(x[-2:])) for x in mat_classification["Node"][z]['ExpirationDate']]
                    Mat['contract_format'][StockType.h300] += [x[2:6] for x in mat_classification["Node"][z]['ExpirationDate']]
        if data['Instruments']["Node"][i]['ENG'] == 'SZSE(O)':
            for mat_classification in data['Instruments']["Node"][i]["Node"][0]["Node"][-4 : ]:
                for z in range(2):
                    QuoteID += mat_classification["Node"][z]['Contracts'] # s300
                    Mat['calendar'][StockType.s300] += [calendar.datetime.date(int(x[0:4]), int(x[4:6]), int(x[-2:])) for x in mat_classification["Node"][z]['ExpirationDate']]
                    Mat['contract_format'][StockType.s300] += [x[2:6] for x in mat_classification["Node"][z]['ExpirationDate']]
        if data['Instruments']["Node"][i]['ENG'] == 'CFFEX(O)':
            for mat_classification in data['Instruments']["Node"][i]["Node"][0]["Node"][-6 : ]:
                for z in range(2):
                    QuoteID += mat_classification["Node"][z]['Contracts'] # gz300
                    Mat['calendar'][StockType.gz300] += [calendar.datetime.date(int(x[0:4]), int(x[4:6]), int(x[-2:])) for x in mat_classification["Node"][z]['ExpirationDate']]
                    Mat['contract_format'][StockType.gz300] += [x[2:6] for x in mat_classification["Node"][z]['ExpirationDate']]


    data = g_QuoteZMQ.QueryAllInstrumentInfo(g_QuoteSession, "Future")
    for i in range(len(data['Instruments']["Node"])):
        if data['Instruments']["Node"][i]['ENG'] == 'CFFEX':
            QuoteID += data['Instruments']["Node"][i]["Node"][2]['Contracts'][1:]
            QuoteID += data['Instruments']["Node"][i]["Node"][3]['Contracts'][1:]
            Mat['calendar'][FutureType.IF] += [calendar.datetime.date(int(x[0:4]), int(x[4:6]), int(x[-2:])) for x in data['Instruments']["Node"][i]["Node"][2]['ExpirationDate'][1:]]
            Mat['contract_format'][FutureType.IF] += [x[2:6] for x in data['Instruments']["Node"][i]["Node"][2]['ExpirationDate'][1:]]
            Mat['calendar'][FutureType.IH] += [calendar.datetime.date(int(x[0:4]), int(x[4:6]), int(x[-2:])) for x in data['Instruments']["Node"][i]["Node"][3]['ExpirationDate'][1:]]
            Mat['contract_format'][FutureType.IH] += [x[2:6] for x in data['Instruments']["Node"][i]["Node"][3]['ExpirationDate'][1:]]


    for sty in [StockType.etf50, StockType.h300, StockType.gz300, StockType.s300, FutureType.IF, FutureType.IH]:
        for format in ['calendar', 'contract_format']:
            Mat[format][sty] = sorted(set(Mat[format][sty]))
            copy = Mat[format][sty].copy()
            Mat[format][sty] = {}
            if sty == StockType.gz300:
                for i, mat in enumerate([Maturity.M1, Maturity.M2, Maturity.M3, Maturity.Q1, Maturity.Q2, Maturity.Q3]):
                    Mat[format][sty][mat] = copy[i]
            elif sty in [StockType.etf50, StockType.h300, StockType.s300, FutureType.IF, FutureType.IH]:
                for i, mat in enumerate([Maturity.M1, Maturity.M2, Maturity.Q1, Maturity.Q2]):
                    Mat[format][sty][mat] = copy[i]

    for sty in [StockType.etf50, StockType.h300, StockType.s300, StockType.gz300]:
        data_opt[sty] = OptData(sty)
        data_opt[sty].subscribe_init(Maturity.M1)
        data_opt[sty].subscribe_init(Maturity.M2)
        data_opt[sty].subscribe_init(Maturity.Q1)
        data_opt[sty].subscribe_init(Maturity.Q2)
        if sty == StockType.gz300:
            data_opt[sty].subscribe_init(Maturity.M3)
            data_opt[sty].subscribe_init(Maturity.Q3)
    for fty in [FutureType.IF, FutureType.IH]:
        data_opt[fty] = FutureData(fty)

    for i in QuoteID + ['TC.S.SSE.510050', 'TC.S.SSE.510300', 'TC.S.SZSE.159919', 'TC.S.SSE.000300']:
        g_QuoteZMQ.subquote(g_QuoteSession,i)


def quote_sub_th(obj,sub_port,filter = ""):
    socket_sub = obj.context.socket(zmq.SUB)
    #socket_sub.RCVTIMEO=7000   #ZMQ超时时间设定
    socket_sub.connect("tcp://127.0.0.1:%s" % sub_port)
    socket_sub.setsockopt_string(zmq.SUBSCRIBE,filter)
    while True:

        message = (socket_sub.recv()[:-1]).decode("utf-8")
        index =  re.search(":",message).span()[1]  # filter
        message = message[index:]
        message = json.loads(message)
        #for message in messages:
        if(message["DataType"]=="REALTIME"):
            server.get_and_cal(message["Quote"])
        elif(message["DataType"] == "PING"):
            g_QuoteZMQ.QuotePong(g_QuoteSession)

    return


def Terminal():
    while True:
        t = time.localtime()
        if t.tm_hour == 15 and t.tm_min > 1:
            g_QuoteZMQKeepAlive.Close()
            server.send('end')


def main():

    global server
    server = Server()

    global g_QuoteZMQKeepAlive
    g_QuoteZMQKeepAlive = KeepAliveHelper(q_data["SubPort"], g_QuoteSession, g_QuoteZMQ)

    #quote
    t = threading.Thread(target = quote_sub_th,args=(g_QuoteZMQ,q_data["SubPort"],))
    t.setDaemon(True)
    t.start()

    t0 = threading.Thread(target = Terminal,args=())
    t0.setDaemon(True)
    t0.start()


if __name__ == '__main__':
    sub_all_options()
    main()
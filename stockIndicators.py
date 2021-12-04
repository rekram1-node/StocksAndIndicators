import yfinance as yf
from googleapiclient.discovery import build
from google.oauth2 import service_account
from statistics import mean
import datetime
import pandas as pd


# This program uses pandas to display the name of each stock listed in a google sheets file that I can update from any device
# it also displays all of the indicators listed below for each stock


# Crossovers (golden cross / death cross)   COMPLETED
# Bollinger Band                            COMPLETED
# MACD                                      COMPLETED
# RSI                                       COMPLETED
# Stochastic Oscillator                     COMPLETED
# Rate of Change                            COMPLETED
# Money Flow Index - (MFI)                  COMPLETED
# Commodity Channel Index (CCI)             COMPLETED




def stockinfo(ticker):
    stock = yf.Ticker(ticker)
    #old  =  stock.history() 

    tod = datetime.datetime.now()



    yearago = datetime.timedelta(days = 365)
    b = tod - yearago
    tod = str(tod).split()
    tod = tod[0]
    b = str(b).split()
    b = b[0]
    old  =  stock.history(start=str(b),  end=str(tod))
    return old


def rsi(info):
    period = 14
    closingPrices = info['Close'][-period-1:]
    changes = []

    for index, close in enumerate(closingPrices):
        if index != 0:
            change = close - closingPrices[index-1]
            changes.append(change)

    positiveChanges = [chng for chng in changes if chng >= 0]
    negativeChanges = [-1 * chng for chng in changes if chng < 0] # positive values are stored but they are negative changes

    AvgPositiveChanges = sum(positiveChanges) / len(positiveChanges)
    AvgNegativeChanges = sum(negativeChanges) / len(negativeChanges)

    RS = AvgPositiveChanges / AvgNegativeChanges
    RSI = 100-(100/(1+RS))

    indicator = 'Uncertain'
    if RSI >= 70:
        indicator = 'Overbought'
    elif RSI <= 30:
        indicator = 'Oversold'



    return "{:.2f}".format(RSI), indicator



def crossover(info):


    data = info

    data1 = data['Close'][-200:]
    Prices = [price for price in data1]
    

    twohundredDayMA = "{:.2f}".format(mean(Prices))
    fiftyDayMA = "{:.2f}".format( mean(Prices[50:len(Prices)]))
    indicator = ''
    if twohundredDayMA > fiftyDayMA:
        indicator = 'Death Cross'
    elif fiftyDayMA > twohundredDayMA:
        indicator = 'Golden Cross'

    return(twohundredDayMA, fiftyDayMA, indicator)  
    

def BollingerBand(info):
    
    # First, calculate a simple moving average. Next, calculate the standard deviation over the 
    # same number of periods as the simple moving average. For the upper band, add the standard deviation 
    # to the moving average. For the lower band, subtract the standard deviation from the moving average.

    Price20 = info['Close'][-20:] # middle band info
    MA = mean(Price20)

    # Upper Band = 20-day SMA + (20-day standard deviation of price x 2) 
    # simple moving average of past 20 days, Middle Band = 20-day simple moving average (SMA)
    # Lower Band = 20-day SMA - (20-day standard deviation of price x 2)

    #  Subtract the moving average from each of the individual data points used in the moving average calculation. 
    # This gives you a list of deviations from the average. Square each deviation and add them all together. 
    # Divide this sum by the number of periods you selected.
    
    deviations = []
    for price in Price20:
        deviations.append((price - MA) ** 2)

    d = sum(deviations) / 20
    #Take the square root of d. This gives you the standard deviation.
    standardDev = d ** 0.5

    upper = MA + standardDev * 2
    lower = MA - standardDev * 2
    middle = MA

    width = upper - lower

    return (width)


def MACD(info):

    
    # Moving Average Convergence Divergence (MACD)
    # Moving average convergence divergence (MACD) is calculated by subtracting the 26-period exponential moving average (EMA) from the 12-period EMA.
    # MACD triggers technical signals when it crosses above (to buy) or below (to sell) its signal line.
    # The speed of crossovers is also taken as a signal of a market is overbought or oversold.
    # MACD helps investors understand whether the bullish or bearish movement in the price is strengthening or weakening.

    twentysixPeriod = info['Close'][-26:] # using the past 26 days of data we take every day and add it to use for a 26 period ema 
    twelvePeriod = info['Close'][-12:]  # using the past 12 days of data we take every day and add it to use for a 12 period ema 

    # K is the weighting factor, the equation is: K = 2 / (n+1) where n is the period
    k26 = 2 / (26 + 1) 
    k12 = 2 / (12 + 1)

    EMA26 = []
    EMA12 = []
   
    # EMA = K * (Current Price - Previous EMA) + Previous EMA


    for index1, price1 in enumerate(twentysixPeriod):
        if index1 == 0: 
            EMA26.append(price1)
        else:
            EMA1 = k26 * (price1 - EMA26[-1]) + EMA26[-1]
            EMA26.append(EMA1)

    for index2, price2 in enumerate(twelvePeriod):
        if index2 == 0: 
            EMA12.append(price2)
        else:
            EMA2 = k12 * (price2 - EMA12[-1]) + EMA12[-1]
            EMA12.append(EMA2)


    # MACD=12-Period EMA − 26-Period EMA
    macd = EMA12[-1] - EMA26[-1]
    #return list(twentysixPeriod), list(twelvePeriod)
    return macd # if values aren't like perfect it is likely due to small rounding error or yfinance not having 100% the best data


def macd(info):
    twelveDays = info['Close'][-13:] # 13 because 13 days ago we will use the sma to calculate ema and use 12 days of emas
    twentySixDays = info['Close'][-27:] # 27 because 27 days ago we will use the sma to calculate ema and use 26 days of emas

    sma12 = [close for close in twelveDays] # not finalized sma
    sma12 = (13 - sum(sma12)) / len(sma12)

    Weightedmultiplier12 = 2/(12+1)
    Weightedmultiplier26 = 2/(26+1)

    sma26 = [close1 for close1 in twentySixDays] # not finalized sma
    sma26 = (27 - sum(sma26) )/ len(sma26)

   # EMA=Price(t)×k+EMA(y)×(1−k)

    EMA12 = []
    EMA26 = []

    for index, price in enumerate(twelveDays):
        if index == 1:
            ema = price * Weightedmultiplier12 + sma12 * (1 - Weightedmultiplier12)
            EMA12.append(ema)
        elif index != 0:
            ema = price * Weightedmultiplier12 + EMA12[index-2] * (1- Weightedmultiplier12)
            EMA12.append(ema)

    for index2, price2 in enumerate(twentySixDays):
        if index2 == 1:
            ema2 = price2 * Weightedmultiplier26 + sma26 * (1 - Weightedmultiplier26)
            EMA26.append(ema2)
        elif index2 != 0:
            ema2 = price2 * Weightedmultiplier26 + EMA26[index2-2] * (1- Weightedmultiplier26)
            EMA26.append(ema2)

    MACDval = EMA12[-1] - EMA26[-1]
    #return EMA12, EMA26
    return MACDval


def StochasticOccilator(info):

    Close = info['Close'][-1]
    Highs = [high for high in info['High'][-14:]]
    Lows = [low for low in info['Low'][-14:]]
    

    Lowest = min(Lows)
    Highest = max(Highs)

    Kprecent = (Close - Lowest) / (Highest - Lowest) * 100
    #     The stochastic oscillator is included in most charting tools and can be easily employed in practice. 
    # The standard time period used is 14 days, though this can be adjusted to meet specific analytical needs. 
    # The stochastic oscillator is calculated by subtracting the low for the period from the current closing price, 
    # dividing by the total range for the period and multiplying by 100. As a hypothetical example, if the 14-day high is $150, 
    # the low is $125 and the current close is $145, then the reading for the current session would be: (145-125) / (150 - 125) * 100, or 80.
    # By comparing the current price to the range over time, the stochastic oscillator reflects the consistency with which the price closes 
    # near its recent high or low. A reading of 80 would indicate that the asset is on the verge of being overbought.
    
    # %K = (Current Close - Lowest Low)/(Highest High - Lowest Low) * 100
    # %D = 3-day SMA of %K

    # Lowest Low = lowest low for the look-back period
    # Highest High = highest high for the look-back period
    # %K is multiplied by 100 to move the decimal point two places
    
    # %K = ((H14−L14) /  (C−L14)) * 100

    # <20 means oversold
    # >80 means overbought
    return Kprecent


def ROC(info):
    
    #     How to Calculate the Price Rate of Change Indicator
    # The main step in calculating the ROC, is picking the "n" value. 
    # Short-term traders may choose a small n value, such as nine. Longer-term investors may 
    # choose a value such as 200. The n value is how many periods ago the current price is being compared to.
    # Smaller values will see the ROC react more quickly to price changes, but that can also mean more false signals. 
    # A larger value means the ROC will react slower, but the signals could be more meaningful when they occur.

    Period = 12
    closing = info['Close'][-1]
    closingPeriodAgo = info['Close'][-Period]

    roc = ((closing - closingPeriodAgo) / closingPeriodAgo) * 100 # as a percent
   # roc2 = ((closing - closingPeriodAgo) / closingPeriodAgo) # not as a percent

    return roc 


def MFI(info):
    0
    #     The Money Flow Index (MFI) is a technical indicator that generates overbought or oversold signals using both prices and volume data.
    # An MFI reading above 80 is considered overbought and an MFI reading below 20 is considered oversold, although levels of 90 and 10 are also used as thresholds.
    # A divergence between the indicator and price is noteworthy. For example, if the indicator is rising while the price is falling or flat, the price could start rising.
  
    p = 14 # number of Periods          <---------- change period accordingly
    closes = info['Close'][-p:]
    highs = info['High'][-p:]
    lows = info['Low'][-p:]
    volumes = info['Volume'][-p:]
   # (High + Low + Close) / 3 = Typical Price
    TypicalPrices = [(high + low +close) /3 for high, low, close in zip(highs, lows, closes)]


    #     Positive Money Flow is calculated by summing the Money Flow of all of the days in the period where Typical Price is higher than the previous period Typical Price.

    # Negative Money Flow is calculated by summing the Money Flow of all of the days in the period where Typical Price is lower than the previous period Typical Price.
    # Typical Price x Volume = Raw Money Flow
    PositiveMoney = []
    NegativeMoney = []
    for index, price in enumerate(TypicalPrices):
        if price > TypicalPrices[index-1] and index != 0:
            RawMoney = price * volumes[index]
            PositiveMoney.append(RawMoney)
        elif index == 0:
            RawMoney = price * volumes[index]
            PositiveMoney.append(RawMoney)
        elif price < TypicalPrices[index-1]:
            RawMoney = price * volumes[index]
            NegativeMoney.append(RawMoney)

    MoneyFlowRatio = sum(PositiveMoney) / sum(NegativeMoney)
   # Money Flow = Typical Price * Volume
    mfi = 100 - 100 / (1 + MoneyFlowRatio)
    return mfi

        # I think this value is somewhat accurate if not I might have slightly misunderstood the equation



    # First, the period's Typical Price is calculated.
    # Typical Price = (High + Low + Close)/3
    # Next, Money Flow (not the Money Flow Index) is calculated by multiplying the period's Typical Price by the volume.
    # Money Flow = Typical Price * Volume
    # If today's Typical Price is greater than yesterday's Typical Price, it is considered Positive Money Flow. If today's price is less, it is considered Negative Money Flow.
    # Positive Money Flow is the sum of the Positive Money over the specified number of periods.
    # Negative Money Flow is the sum of the Negative Money over the specified number of periods.
    # The Money Ratio is then calculated by dividing the Positive Money Flow by the Negative Money Flow.
    # Money Ratio = Positive Money Flow / Negative Money Flow
    # Finally, the Money Flow Index is calculated using the Money Ratio



def CCI(info):
    # Determine how many periods your CCI will analyze. Twenty is commonly used. Fewer periods result in a 
    # more volatile indicator, while more periods will make it smoother. For this calculation, we will assume 20 periods.
    # Adjust the calculation if using a different number.

    p = 20 # number of Periods
    closes = info['Close'][-p:]
    highs = info['High'][-p:]
    lows = info['Low'][-p:]

    TypicalPrices = []

    for c, h, l in zip(closes, highs, lows):
        typicalP = (h + l + c) / 3
        TypicalPrices.append(typicalP)

    MA = sum(TypicalPrices) / p
    TypicalPrice = TypicalPrices[-1]

    #     There are four steps to calculating the Mean Deviation: 
    # First, subtract the most recent 20-period average of the typical price from each period's typical price. 
    # Second, take the absolute values of these numbers. 
    # Third, sum the absolute values. 
    # Fourth, divide by the total number of periods (20). 

    meanDev = sum([abs(i - MA) for i in TypicalPrices]) / p

    cci = (TypicalPrice - MA) / (0.015 * meanDev)

    #cci = (TypicalPrices - MA) / (0.015 * MeanDeviation)
    # When the CCI moves from negative or near-zero territory to above 100, 
    # that may indicate the price is starting a new uptrend. Once this occurs, 
    # traders can watch for a pullback in price followed by a rally in both price 
    # and the CCI to signal a buying opportunity.



    return cci




def Names_Tickers (range):
    RANGE = range # need a way to figure out when the range gets bigger  # 'A1:B3'

    SCOPES = ['https://www.googleapis.com/auth/spreadsheets.readonly'] # delete the .readonly if you want to be able to write to it as well

    SERVICE_ACCOUNT_FILE = 'stockkeys.json'

    creds = None
    creds = service_account.Credentials.from_service_account_file(
        SERVICE_ACCOUNT_FILE, scopes = SCOPES
    )


    # The ID and range of a sample spreadsheet.
    SAMPLE_SPREADSHEET_ID = '1Gt3NxokBIId3W3zD_0SGDLAW3KN86KK6E5NG-YUcl30'


    service = build('sheets', 'v4', credentials=creds)

    # Call the Sheets API
    sheet = service.spreadsheets()
    result = sheet.values().get(spreadsheetId=SAMPLE_SPREADSHEET_ID,
                                range=RANGE).execute()    # range needs to be updated in accordance to how many stocks are added!

    values = result.get('values', [])

    return (values)









a = Names_Tickers('d6') # d6 lists the range of columns to read up to on the google sheets
Range = a[0][0]

dataframe = {   # dictionary to use for pandas
    'Price': [],
    'RSI': [],
    'Cross Over': [],
    'Bollinger Band Width':[],
    'MACD':[],
    'Stochastic Occilator':[],
    'ROC':[],
    'MFI':[],
    'CCI':[]
}

names = []

for nameTick in Names_Tickers(Range):
    names.append(nameTick[0])
    Symbol = nameTick[1]

    StockData = stockinfo(Symbol)
    dataframe['Price'].append(StockData['Close'][-1])
    dataframe['RSI'].append(rsi(StockData))
    dataframe['Cross Over'].append(crossover(StockData))
    dataframe['Bollinger Band Width'].append(BollingerBand(StockData))
    dataframe['MACD'].append(MACD(StockData))
    dataframe['Stochastic Occilator'].append(StochasticOccilator(StockData))
    dataframe['ROC'].append(ROC(StockData))
    dataframe['MFI'].append(MFI(StockData))
    dataframe['CCI'].append(CCI(StockData))





df = pd.DataFrame(dataframe, names)

print(df) 






####################### DONE ##############################
# MFI seems right on the money
# Cross over is correct
# CCI seems to be perfect
# ROC is fairly close
# Stochastic Occilator has the exact equation so it preforms as intended though it may through a slightly different value than some other platforms
# bollinger band width has the exact equation, the prices seem to be accurate but it differs from some platforms

###################### To Fix #############################

# MACD   might be messed up

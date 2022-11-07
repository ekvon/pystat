# Created:2018-12-13
# Modified:2021-07-29
# Part Of:libstat Package
# Description:Realization of histogram concept.
import math
from . import axis
from .axis import IntervalAxis
class Interval:
    def __init__(self):
        'Empty histogram'
        self._axis=IntervalAxis()
        # storage
        self._freq=[]

    def shape(self):
        'Return number of intervals in the series'
        return tuple(self._axis.size())

    def regular(self,min_,max_,num_parts):
        'Split on intervals uniformly'
        # parameters are checked in axis-method
        self._axis.regular(min_,max_,num_parts)
        self._freq=[0 for i in range(num_parts)]

    def sturges(self,min_,max_,volume):
        'Split on intervals using Sturges rule'
        self._axis.sturges(min_,max_,volume)
        size=self._axis.size()
        if size<=0:
            return
        self._freq=[0 for i in range(size)]

    def variable(self,min_,max_,num_parts,fn):
        'Split on intervals with the help of callable object'
        self._axis.variable(min_,max_,num_parts,fn)
        size=self._axis.size()
        if size<=0:
            return
        self._freq=[0 for i in range(size)]

    # @param x (iterable) - loaded values
    def load(self,x):
        'Initialization with the help of observation series'
        N=len(x)
        if N==0:
            return
        for i in range(N):
            idx=self._axis.index(x[i])
            if idx<0:
                # unexpected
                continue
            # increment frequency
            self._freq[idx]+=1

    def index(self,x):
        'Retrieve cell in which specified attribute value is located'
        return self._axis.index(x)

    # @param idx (iterable) - index of multidimension interval
    def freq(self,idx):
        'Retrieve frequency of specifiedinterval'
        idx=int(idx)
        size=len(self._freq)
        if idx<0 or size<=idx:
            return -1
        return self._freq[idx]

    def range(self,idx):
        'Return the range of specified interval'
        return self._axis.range(idx)

    def bounds(self):
        'Return array of numbers which are bounds of intervals'
        return self._axis._data
    
    def size(self):
        'Retrieve number of intervals in decomposition'
        return self._axis.size()

    def volume(self):
        'Return frequency of whole series'
        return sum(self._freq)

    def freq(self):
        'Return array of intervals frequencies'
        return tuple(self._freq)

    def median(self):
        'Return index of median interval and median estimation'
        # The following formula is used to get estimation of median:
        #   Me=a+h*(0.5-cum_freq)/med_freq, where
        #   a - the begining of median interval
        #   h - the range of median interval
        #   cum_freq - cummulative frequency (related) before median interval
        #   med_freq - cummulative frequency (related) of median interval
        # In detail see formula 1.16 at p.30.
        # Total number of observations
        N=self.volume()
        median=divmod(N,2)[0]
        # cummulative frequency
        cummulative=0
        i=0
        size=self.size()
        # definition of interval in which median value is placed
        while(i<size):
            freq=self._freq[i]
            if(median<=(cummulative+freq)):
                break
            else:
                cummulative+=freq
                # next (i)
                i+=1
        r=self._axis.range(i)
        # the length of median interval
        h=r[1]-r[0]
        Wc=cummulative/N
        Wm=freq/N
        result=r[0]+h*(0.5*Wc)/Wm
        return i,result

    def mode(self):
        'Return index of mode interval and mode value estimation'
        # Mode interval is interval with the highest frequency.
        # In general, for mode value estimation the following formula is used:
        #   mode=a+h*(Wm-Wp)/(2*Wm-Wp-Wn)
        #   a - the begining of mode interval
        #   h - mode interval range
        #   Wm - relative frequency of mode interval
        #   Wp - relative frequency of previuos interval
        #   Wn - relative frequency of next interval
        # in detail, see formula 1.18 at p.31.
        N=self.volume()
        max_freq=max(self._freq)
        idx=self._freq.index(max_freq)
        # mode interval size
        r=self._axis.range(idx)
        h=r[1]-r[0]
        # estimation
        if idx==0:
            Wm=self._freq[idx]/N
            Wn=self._freq[1]/N
            return idx,r[0]+h*(Wm-Wn)/(2*Wm-2*Wn)
        elif idx==(N-1):
            Wm=self._freq[idx]/N
            Wpn=self._freq[idx-1]/N
            return idx,r[0]+h*(Wm-Wp)/(2*Wm-2*Wp)
        else:
            Wm=self._freq[idx]/N
            Wp=self._freq[idx-1]/N
            Wn=self._freq[idx+1]/N
            return idx,r[0]+h*(Wm-Wp)/(2*Wm-Wp-Wn)
        
    def F(self,i,absolute=True):
        'Empirical distribution function.'
        # Return cummulative frequency for specified interval. That's the
        # the sum of frequencies (absolute or relative)of previous intervals.
        # If index 'i' is '0' then '0' is returned.
        # If index 'i' equals to the size of series (number of intervals) then
        # the frequency of whole series is returned.
        # For non correct value of index the None is returned.
        size=len(self.m_Intervals)
        if(type(i)!=int or i<0 or size<i):
            return None
        result=0
        for k in range(0,i):
            result+=self.m_Intervals[k][2]
        if(absolute==True):
            return result
        frequency=self.Freq()
        return result/frequency

    def p(self,x):
        ''
        index=self.GetByValue(x)
        if index==None:
            return 0
        freq=self.Freq()
        return self.FreqOfInterval(index)/self.Freq()
        
def rejection(values,threshold=0.9,repeat=True):
    ''
    # @param values - list of real or integers.
    # Description:
    # Interval series is created. If relative frequency of some interval is
    # more than threshold then all values outside this interval are removed
    # from initial list. Recursive call is not realized yet.
    s=Interval(values)
    count=s.GetCount()
    freq=s.Freq()
    for i in range(0,len(freq)):
        if threshold<freq[i]/count:
            # (inf,sup,freq)
            t=s.GetByIndex(i)
            k=0
            # modification of initial values
            size=len(values)
            while k<size:
                if values[k]<t[0] or t[1]<values[k]:
                    # index is unchanged
                    values.pop(k)
                else:
                    k+=1
            return True,values
        else:
            pass
    # initial values are unchanged
    return False,values

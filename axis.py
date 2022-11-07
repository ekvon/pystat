# Created:2021-07-31
# Type:Python Module
# Part Of:libstat Package
import math
class IntervalAxis:
    def __init__(self):
        ''
        # empty
        self._data=()

    def __call__(self,x):
        ''
        return self.index(x)

    # Split interval uniformly on specified number of parts
    # @param min_ - the lower bound of interval
    # @param max_ - the upper bound of interval
    # @param num_parts - number of parts
    def regular(self,min_,max_,num_parts):
        ''
        fMin=float(min_)
        fMax=float(max_)
        N=int(num_parts)
        if fMax<=fMin:
            raise ValueError('Invalid interval bounds')
        if N<=1:
            raise ValueError('Invalid number of parts')
        h=(fMax-fMin)/N
        data=[fMin+i*h for i in range(N+1)]
        self._data=tuple(data)

    # Split interval using Sturges' rule.
    # @param volume - the length of observation series.
    def sturges(self,min_,max_,volume):
        'Initialization with the help of Sturges rule'
        fMin=float(min_)
        fMax=float(max_)
        vol=int(volume)
        if fMax<=fMin:
            raise ValueError('Invalid interval bounds')
        if volume<=1:
            raise ValueError('Invalid observation series length')
        data=[]
        h=(fMax-fMin)/(1+3.322*math.log10(vol))
        a=fMin-h/2
        i=0
        while a+h*i<fMax:
            data.append(a+h*i)
            # next (i)
            i+=1
        # the last value
        data.append(a+h*i)
        self._data=tuple(data)

    def variable(self,min_,num_parts,fn):
        'Split interval using callable object'
        fMin=float(min_)
        N=int(num_parts)+1
        if N<=1:
            raise ValueError('Invalid number of parts')
        data=[fn(0)]
        for n in range(1,N):
            val=fn(n)
            if val<=data[n-1]:
                raise ValueError('The sequence is breaked')
            data.append(val)
        self._data=tuple(data)

    def size(self):
        ''
        return len(self._data)-1

    # Return index of cell where 'x' is located
    # Params:
    # @x (float) - value for which cell is searched.
    # Remarks:
    # Linear search is used.
    def index(self,x):
        ''
        N=len(self._data)-1
        i=-1
        while i<N:
            if x<self._data[i+1]:
                return i
            # next (i)
            i+=1
        # not found
        return -1

    def range(self,idx):
        ''
        idx=int(idx)
        num_intervals=len(self._data)-1
        if idx<0 or num_intervals<=idx:
            raise ValueError('invalid index value: {}'.format(idx))
        return (self._data[idx],self._data[idx+1])


class vector2:

    def __init__(self, x0, y0, x1, y1):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1

    def __init__(self, x1, y1, x0 = 0, y0 = 0):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1

    @property
    def median_point(self):
        return  vector2((self.x0+self.x1)/2, (self.y0+self.y1)/2)

    @property
    def magnitude(self):
        return  ((self.x1-self.x0)**2+(self.y1-self.y0)**2)**0.5
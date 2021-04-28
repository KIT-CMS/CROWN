class Quantity:
    def __init__(self, name):
        self.name = name
        self.shifts = []
        self.children = []

    def get_leaf(self, shift):
        # add sth here that creates name of shifted quantity
        if shift in self.shifts:
            return self.name + shift
        return self.name

    def shift(self, name):
        if not name in self.shifts:
            self.shifts.append(name)
            for c in self.children:
                c.shift(name)


pt_1 = Quantity("pt_1")
pt_2 = Quantity("pt_2")
m_vis = Quantity("m_vis")
mt_1 = Quantity("mt_1")

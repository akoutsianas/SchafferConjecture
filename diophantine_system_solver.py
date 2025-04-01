class DiophantineSystem:

    def __init__(self, f, d1, d2, n_greater_than=None):
        self.f = f
        self.d1 = d1
        self.d2 = d2
        self.k = f.degree()
        self._ais = self._compute_ai()
        if self._ais[0] * self.d1 * self.d2**self.k < 0:
            raise ValueError(f"The condition a0*d1*d2^k > 0 is not satisfied!")
        self.k0 = self._compute_k0()
        self._D = (self._ais[0] * self.d2**self.k) / self.d1
        self._theta0 = self._compute_theta0()
        self.m = self._compute_m()
        self.C0 = None
        self.B2 = None
        self.alpha = None

    def _compute_ai(self):
        ais = [self.f.coefficient(self.k - i) for i in range(self.k + 1)]
        return ais

    def _compute_k0(self):
        k0 = min([i for i in range(1, self.k + 1) if self._ais[i] != 0])
        return k0

    def _compute_theta0(self):
        if self._D > 1:
            return 1
        else:
            # TODO we may choose a bigger choice of theta0 for bigger n instead of n=2 and give a lower bound of n
            return sqrt(self._D)

    def bound_n(self):
        self.C0 = self._compute_C0()

    def _compute_m(self):
        if self._D < 1:
            return 1
        else:
            return round(self._D)

    def _compute_C0(self):
        pass

class DiophantineSystem:

    def __init__(self, f, d1, d2, prec=20, n_greater_than=None):
        self._RR = RealField(prec=prec)
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
        self.C1 = self._compute_C1()
        self.C2 = self._compute_C2()
        self.B2 = None
        self.alpha = None
        self._a1 = None

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

    def _compute_m(self):
        if self._D < 1:
            return 1
        else:
            return round(self._D)

    def _compute_C1(self):
        s2 = self.d2.denominator()
        return 1/(s2**(self.k - self.k0))

    def _compute_C2(self):
        return 1 + max(1, self._theta0)

    def bound_n(self):
        self.C0, self.B2, self.alpha = self._compute_suitable_C0_B2_alpha()

    def _compute_suitable_C0_B2_alpha(self, p=0.1):
        self._a1 = abs(self._D - self.m)
        cis = [self.d2**(self.k - i) * self._ais[i] / self.d1 for i in range(self.k0, self.k+1)]
        print(f"cis: {cis}")
        B2 = 2
        C0 = sum([abs(ci) / B2**i for i, ci in enumerate(cis)])
        alpha = self._compute_alpha(C0, B2)
        print(f"Initial, C0: {C0}, alpha: {alpha}")
        self._lower_bound_for_B2_is_satisfied(B2, C0, alpha)
        while ((abs(alpha - self._a1) / self._a1) >= p) or (not self._lower_bound_for_B2_is_satisfied(B2, C0, alpha)):
            B2 += 1
            C0 = sum([abs(ci) / B2 ** i for i, ci in enumerate(cis)])
            alpha = self._compute_alpha(C0, B2)
        print(f"C0: {C0}, B2: {B2}")
        return C0, B2, alpha

    def _lower_bound_for_B2_is_satisfied(self, B2, C0, alpha):
        if alpha >= 1:
            return False
        rt1 = (self._RR(C0/self._theta0)).nth_root(self.k0)
        rt2 = (self._RR(C0/(1 - alpha))).nth_root(self.k0)
        return B2 > max(C0 + 2, rt1, rt2)

    def _compute_alpha(self, C0, B2):
        """
        We compute alpha either from its definition \alpha^\prime + C/B_2^{k_0} or we choose alpha to be equal
        alpha^\prime + (1 - p)\alpha^\\prime for some given p < 1.
        """
        return self._a1 + C0 / B2**self.k0

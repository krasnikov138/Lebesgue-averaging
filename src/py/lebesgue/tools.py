import pandas as pd
from lebesgue.bindings import PyEndfFile, PyInterpolationTable


class InterpolationTable:
    def __init__(self, cdata, namex, namey):
        if not isinstance(cdata, PyInterpolationTable):
            raise ValueError("cdata obj must be PyInterpolationTable")
        self._cdata = cdata
        self.namex = namex
        self.namey = namey

    @property
    def interpolation_schema(self):
        return pd.DataFrame({
            'boundary': self._cdata.boundaries,
            'interpolation_type': self._cdata.interpolation_types
        })

    @property
    def data(self):
        return pd.DataFrame({
            self.namex: self._cdata.xs,
            self.namey: self._cdata.ys,
        })

    def __len__(self):
        return self._cdata.xs.shape[0]

    def __repr__(self):
        return "Interpolation Schema\n" + self.interpolation_schema.__repr__() +\
               "\nInterpolation Table\n" + self.data.__repr__()


class EndfData:
    def __init__(self, payload):
        for key in payload.keys():
            setattr(self, key, payload[key])


class CrossSectionData(EndfData):
    def __init__(self, payload):
        super().__init__(payload)
        self.cs = InterpolationTable(self.cs, "energy", "cs")

    def __repr__(self):
        return f"CrossSectionData(ZA={self.za}, AWR={self.awr:.2f}, QM={self.qm:.2f}, QI={self.qi:.2f}, LR={self.lr})"


class ContinuumDistribution(EndfData):
    def __repr__(self):
        string_repr = f"ContinuumDistribution(LANG={self.angular_repr}, " +\
                      f"points={len(self.primary_energies)})"
        return string_repr


class ProductSubsection(EndfData):
    def __init__(self, payload):
        super().__init__(payload)
        self.energy_yield = InterpolationTable(self.energy_yield, "energy", "yield")
        self.distr = ContinuumDistribution(self.distr)

    def __repr__(self):
        attrs = ['ZAP', 'atomic_weight_ratio', 'product_modifier_flag', 'law']
        names = ['ZAP', 'AWR', 'LIP', 'LAW']

        output = [f'{name}={getattr(self, attr):.2f}' for name, attr in zip(names, attrs)]
        return f"ProductSubsection({', '.join(output)})"

    def describe(self):
        string_repr = f"Product charge/mass parameters (ZAP) = {self.ZAP}\n"
        string_repr += f"Atomic weight ratio (AWR) = {self.atomic_weight_ratio}\n"
        string_repr += f"Product modifier flag (LIP) = {self.product_modifier_flag}\n"
        string_repr += f"Distribution law (LAW) = {self.law}"
        print(string_repr)


class EnergyAngleData(EndfData):
    def __init__(self, payload):
        super().__init__(payload)
        self.subsections = [ProductSubsection(section) for section in self.subsections]

    def __repr__(self):
        attrs = ['ZA', 'atomic_weight_ratio', 'yield_multiplicity_flag', 'reference_frame']
        names = ['ZA', 'AWR', 'JP', 'LCT']

        output = [f'{name}={getattr(self, attr):.2f}' for name, attr in zip(names, attrs)]
        return f"EnergyAngleData({', '.join(output)}, subsections={len(self.subsections)})"

    def describe(self):
        string_repr = f"EnergyAngleData object contains {len(self.subsections)} subsections\n"
        string_repr += f"Substance charge and mass parameters (ZA) = {self.ZA}\n"
        string_repr += f"Atomic weight ratio (AWR) = {self.atomic_weight_ratio}\n"
        string_repr += f"Promt fission neutrons and photons (JP) = {self.yield_multiplicity_flag}\n"
        string_repr += f"Reference system for energy and angle (LCT) = {self.reference_frame}"
        print(string_repr)


class EndfFile(PyEndfFile):

    @property
    def table_of_content(self):
        return pd.DataFrame(super().table_of_content, columns=['mf', 'mt'])

    def get_cross_section(self, mt):
        return CrossSectionData(super().get_cross_section(mt))

    def get_energy_angle(self):
        return EnergyAngleData(super().get_energy_angle())

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

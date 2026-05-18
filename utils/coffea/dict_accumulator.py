from coffea import processor


class DictAccumulator(processor.AccumulatorABC):
    def __init__(self, value={}):
        self.value = value

    def identity(self):
        return DictAccumulator({})
    
    def add(self, other):
        if len(other.value.keys()) == 0:
            return
        elif len(self.value.keys()) == 0:
            self.value = other.value
        else:
            if set(other.value.keys()) != set(self.value.keys()):
                missing_in_self = set(other.value.keys()) - set(self.value.keys())
                missing_in_other = set(self.value.keys()) - set(other.value.keys())
                raise ValueError(
                    f"DictAccumulator keys mismatch!\n"
                    f"Keys in other but not in self: {missing_in_self}\n"
                    f"Keys in self but not in other: {missing_in_other}\n"
                    f"Self keys: {set(self.value.keys())}\n"
                    f"Other keys: {set(other.value.keys())}"
                )
            for key in other.value.keys():
                self.value[key] += other.value[key]


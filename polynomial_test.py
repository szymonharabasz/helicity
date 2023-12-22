import torch
import math


class Heliity1d(torch.nn.Module):
    def __init__(self):
        """
        In the constructor we instantiate four parameters and assign them as
        member parameters.
        """
        super().__init__()
        self.lambda_theta = torch.nn.Parameter(torch.randn(()))

    def forward(self, x):
        """
        In the forward function we accept a Tensor of input data and we must return
        a Tensor of output data. We can use Modules defined in the constructor as
        well as arbitrary operators on Tensors.
        """
        return torch.vstack(
            (1.0 + self.lambda_theta * x[0] ** 2) * x[1],
            (1.0 + self.lambda_theta * x[0] ** 2) * x[2])
    def string(self):
        """
        Just like any class in Python, you can also define custom method on PyTorch modules
        """
        return f'y = (1 + {self.lambda_theta.item()} x^2'

# Construct our model by instantiating the class defined above
model = Heliity1d()
for param in model.parameters():
    print(param)

import torch
import torch.nn as nn
import torch.nn.functional as F
from einops import rearrange
from pyannote.audio import Model
from pyannote.audio.models.segmentation import PyanNet
from pyannote.audio.models.segmentation.debug import SimpleSegmentationModel

SNR_MIN = -15
SNR_MAX = 80
C50_MIN = -10
C50_MAX = 60


class ParametricSigmoid(nn.Module):
    def __init__(self, alpha: float, beta: float) -> None:
        super().__init__()
        self.alpha = alpha
        self.beta = beta

    def forward(self, x: torch.Tensor):
        return (self.beta - self.alpha) * F.sigmoid(x) + self.alpha


class CustomClassifier(nn.Module):
    def __init__(self, in_features, out_features: int) -> None:
        super().__init__()
        # Use direct attributes instead of ModuleDict to avoid dictionary overhead in forward pass
        # This provides minor but measurable speedup since forward is called millions of times
        self.vad_head = nn.Linear(in_features, out_features)
        self.snr_head = nn.Linear(in_features, 1)
        self.c50_head = nn.Linear(in_features, 1)

    def forward(self, x: torch.Tensor):
        # Direct calls avoid Python loop overhead - compiler can optimize better
        return {
            'vad': self.vad_head(x),
            'snr': self.snr_head(x),
            'c50': self.c50_head(x)
        }


class CustomActivation(nn.Module):
    # Try something else for snr and c50
    # TODO : print and save alpha and beta values of ParametricSigmoid
    def __init__(self) -> None:
        super().__init__()
        # Use direct attributes to avoid loop overhead in forward pass
        self.vad_act = nn.Sigmoid()
        self.snr_act = ParametricSigmoid(SNR_MIN, SNR_MAX)
        self.c50_act = ParametricSigmoid(C50_MIN, C50_MAX)

    def forward(self, x: torch.Tensor):
        # Direct tensor operations - no list/stack, just concat
        # This avoids list allocation and intermediate stacking operations
        vad_out = self.vad_act(x['vad'])
        snr_out = self.snr_act(x['snr'])
        c50_out = self.c50_act(x['c50'])

        # Use torch.cat instead of stack+rearrange for better performance
        return torch.cat([vad_out, snr_out, c50_out], dim=-1)


class RegressiveSegmentationModelMixin(Model):
    def build(self):
        """
        Debug architecture
        """
        nb_classif = len(set(self.specifications.classes) - set(['snr', 'c50']))
        self.classifier = CustomClassifier(32 * 2, nb_classif)
        self.activation = CustomActivation()


class CustomSimpleSegmentationModel(RegressiveSegmentationModelMixin, SimpleSegmentationModel):
    pass


class CustomPyanNetModel(RegressiveSegmentationModelMixin, PyanNet):
    def build(self):
        if self.hparams.linear["num_layers"] > 0:
            in_features = self.hparams.linear["hidden_size"]
        else:
            in_features = self.hparams.lstm["hidden_size"] * (
                2 if self.hparams.lstm["bidirectional"] else 1
            )
        nb_classif = len(set(self.specifications.classes) - set(['snr', 'c50']))
        self.classifier = CustomClassifier(in_features, nb_classif)
        self.activation = CustomActivation()

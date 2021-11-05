import logging
from typing import List, Union

from code_generation.producer import (
    CollectProducersOutput,
    Producer,
    ProducerGroup,
    TProducerInput,
    TProducerStore,
)
from code_generation.quantity import QuantitiesStore

log = logging.getLogger(__name__)


class ProducerRule:
    def __init__(
        self,
        producers: TProducerInput,
        samples: Union[str, List[str]],
        scopes: Union[str, List[str]] = "global",
        invert: bool = False,
        update_output: bool = True,
    ):
        if isinstance(producers, ProducerGroup) or isinstance(producers, Producer):
            producers = [producers]
        self.producers = producers
        if isinstance(samples, str):
            samples = [samples]
        self.samples = samples
        if isinstance(scopes, str):
            scopes = [scopes]
        self.scopes = scopes
        self.invert = invert
        self.update_output = update_output
        self.global_scope = "global"

    def set_scopes(self, scopes: List[str]) -> None:
        if isinstance(scopes, str):
            scopes = [scopes]
        self.scopes = scopes

    def set_global_scope(self, global_scope: str) -> None:
        self.global_scope = global_scope

    # Evaluate whether modification should be applied depending on sample and inversion flag
    def is_applicable(self, sample: str) -> bool:
        applicable = sample in self.samples
        if self.invert:
            applicable = not applicable
        return applicable

    # Placeholder for the actual operation on a list. To be overwritten by inheriting classes
    def update_producers(self, producers_to_be_updated: TProducerStore) -> None:
        log.error("Operation not implemented for ProducerRule base class!")
        pass

    def update_outputs(self, outputs_to_be_updated: QuantitiesStore) -> None:
        log.error("Operation not implemented for ProducerRule base class!")
        pass

    def apply(
        self,
        sample: str,
        producers_to_be_updated: TProducerStore,
        outputs_to_be_updated: QuantitiesStore,
    ) -> None:
        if self.is_applicable(sample):
            log.info("For sample {}, applying >> {} ".format(sample, self))
            self.update_producers(producers_to_be_updated)
            self.update_outputs(outputs_to_be_updated)


# Modifier class that can remove producers from lists
class RemoveProducer(ProducerRule):
    def __str__(self) -> str:
        return "ProducerRule - remove {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def __repr__(self) -> str:
        return "ProducerRule - remove {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def update_producers(self, producers_to_be_updated: TProducerStore):
        for scope in self.scopes:
            for producer in self.producers:
                if producer in producers_to_be_updated[scope]:
                    log.debug(
                        "RemoveProducer: Removing {} from producer in scope {}".format(
                            producer, scope
                        )
                    )
                    producers_to_be_updated[scope].remove(producer)

    def update_outputs(self, outputs_to_be_updated: QuantitiesStore):
        if self.update_output:
            outputs: QuantitiesStore = {}
            # if the producer is in the global scope, we add the output to all running scopes
            if self.scopes == [self.global_scope]:
                log.debug(
                    "Updating outputs for producer in global scope --> adding output to all scopes"
                )
                scopes = outputs_to_be_updated.keys()
                for scope in scopes:
                    outputs[scope] = CollectProducersOutput(
                        self.producers, self.global_scope
                    )
            else:
                scopes = self.scopes
                for scope in scopes:
                    outputs[scope] = CollectProducersOutput(self.producers, scope)
            for scope in scopes:
                for output in outputs[scope]:
                    if output in outputs_to_be_updated[scope]:
                        log.debug(
                            "RemoveProducer: Removing {} from outputs in scope {}".format(
                                output, scope
                            )
                        )
                        outputs_to_be_updated[scope].remove(output)


# Modifier class that can append producers to lists
class AppendProducer(ProducerRule):
    def __str__(self) -> str:
        return "ProducerRule - add {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def __repr__(self) -> str:
        return "ProducerRule - add {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def update_producers(self, producers_to_be_updated: TProducerStore):
        for scope in self.scopes:
            for producer in self.producers:
                log.debug(
                    "AppendProducer: Adding {} to producers in scope {}".format(
                        producer, scope
                    )
                )
                producers_to_be_updated[scope].append(producer)

    def update_outputs(self, outputs_to_be_updated: QuantitiesStore):
        if self.update_output:
            outputs: QuantitiesStore = {}
            # if the producer is in the global scope, we add the output to all running scopes
            if self.scopes == [self.global_scope]:
                log.debug(
                    "Updating outputs for producer in global scope --> adding output to all scopes"
                )
                scopes = outputs_to_be_updated.keys()
                for scope in scopes:
                    outputs[scope] = CollectProducersOutput(
                        self.producers, self.global_scope
                    )
            else:
                scopes = self.scopes
                for scope in scopes:
                    outputs[scope] = CollectProducersOutput(self.producers, scope)
            for scope in scopes:
                for output in outputs[scope]:
                    log.debug(
                        "AppendProducer: Adding {} to outputs in scope {}".format(
                            output, scope
                        )
                    )
                    outputs_to_be_updated[scope].add(output)

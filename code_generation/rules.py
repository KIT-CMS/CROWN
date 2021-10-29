import logging
from config.utility import CollectProducersOutput


log = logging.getLogger(__name__)


class ProducerRule:
    def __init__(
        self, producers, samples, scopes="global", invert=False, update_output=True
    ):
        self.producers = producers
        self.samples = samples
        self.scopes = scopes
        self.invert = invert
        self.update_output = update_output
        self.global_scope = "global"

    def set_scopes(self, scopes):
        self.scopes = scopes

    def set_global_scope(self, scope):
        self.global_scope = scope

    # Evaluate whether modification should be applied depending on sample and inversion flag
    def applicable(self, sample):
        applicable = sample in self.samples
        if self.invert:
            applicable = not applicable
        return applicable

    # Placeholder for the actual operation on a list. To be overwritten by inheriting classes
    def update_producers(self, producer, producers_to_be_updated):
        log.error("Operation not implemented for ProducerRule base class!")
        pass

    def update_outputs(self, output, outputs_to_be_updated):
        log.error("Operation not implemented for ProducerRule base class!")
        pass

    def apply(self, sample, producers_to_be_updated, outputs_to_be_updated):
        if self.applicable(sample):
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

    def update_producers(self, producers_to_be_updated):
        for scope in self.scopes:
            for producer in self.producers:
                if producer in producers_to_be_updated[scope]:
                    log.debug(
                        "RemoveProducer: Removing {} from producer in scope {}".format(
                            producer, scope
                        )
                    )
                    producers_to_be_updated[scope].remove(producer)

    def update_outputs(self, outputs_to_be_updated):
        if self.update_output:
            outputs = {}
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

    def update_producers(self, producers_to_be_updated):
        for scope in self.scopes:
            for producer in self.producers:
                log.debug(
                    "AppendProducer: Adding {} to producers in scope {}".format(
                        producer, scope
                    )
                )
                producers_to_be_updated[scope].add(producer)

    def update_outputs(self, outputs_to_be_updated):
        if self.update_output:
            outputs = {}
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

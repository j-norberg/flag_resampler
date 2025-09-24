#ifndef STREAMER_H
#define STREAMER_H

class ISampleProducer;

ISampleProducer* streamer_factory(ISampleProducer* input, int sr_to);

#endif


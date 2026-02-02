#ifndef STREAMER_H
#define STREAMER_H

class ISampleProducer;

// rule of thumb
// 100% is as high as this algo can go
// 50% is ok if destination is less than 64-bit doubles
ISampleProducer* streamer_factory(ISampleProducer* input, int sr_to, int quality_percentage);

#endif


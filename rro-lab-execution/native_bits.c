/* Each node is a literal ML-I primitive location, with an actual next link. */
#include <stdint.h>
#include <stddef.h>
typedef struct Wire { struct Wire *next; uint16_t byte; uint8_t mask, action, grant; } Wire;
typedef struct Init {
    const uint8_t *source; uint8_t *target; Wire *cursor;
    uint8_t latch, done, grant;
    uint64_t counts[4];
} Init;
uint64_t initialize(Init *s, uint64_t limit) {
    uint64_t n=0;
    while(n<limit && !s->done) {
        Wire *p=s->cursor;
        switch(p->action) {
        case 0: s->target[p->byte] &= (uint8_t)~p->mask; break;
        case 1: s->latch=!!(s->source[p->byte]&p->mask); break;
        case 2:
            if(s->latch) s->target[p->byte] |= p->mask;
            else s->target[p->byte] &= (uint8_t)~p->mask;
            break;
        case 3: s->done=1; break;
        default: return UINT64_MAX;
        }
        s->counts[p->action]++; s->cursor=p->next;
        n++;
    }
    return n;
}

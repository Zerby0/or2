#ifndef LIST_H
#define LIST_H

#define _LIST_DEF(name, type) \
	typedef struct { \
		type* buf; \
		int len; \
		int capacity; \
	} List_ ## name; \
	void list_ ## name ## _init(List_ ## name* l); \
	void list_ ## name ## _push(List_ ## name* l, type val); \
	void list_ ## name ## _pop_front(List_ ## name* l);

_LIST_DEF(d, double)

#endif /* LIST_H */

.globl atomic_add
atomic_add:
    lock xaddq %rsi, (%rdi)
	ret

.globl atomic_get
atomic_get:
    movq (%rdi), %rax
	ret

.globl atomic_set
atomic_set:
    movq %rsi, (%rdi)
	ret
	
.globl atomic_inc
atomic_inc:
    lock incq (%rdi)
    ret


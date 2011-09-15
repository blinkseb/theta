.text
.align 4
.globl log2_dot
.globl template_nllikelihood
log2_dot:
    fldz
    testl %edx, %edx
    je .exit
    subl $1, %edx
    xorl %eax, %eax
    leaq 8(,%rdx,8), %rdx
    .p2align 4
.loopstart:
    fldl (%rsi,%rax)
    fldl (%rdi,%rax)
    fyl2x
    addq $8, %rax
    cmpq %rdx, %rax
    faddp %st, %st(1)
    jne .loopstart
.exit:
    fstpl -8(%rsp)
    movsd -8(%rsp), %xmm0
    ret

/* note: according to the System V AMD64 ABI, we only have to preserve %rbx, %rsp, %rbp, %r12-%r15, so we just do not use these ... */
.p2align 4
template_nllikelihood:  /* data = %rdi, pred = %rsi; n = %edx */
    fldz /* use st(0) to save    sum_i  data[i] * log(pred[i])   */
    xorpd %xmm2, %xmm2 /* use xmm2 to save     sum_i pred[i] */
    testl %edx, %edx
    je .tl_exit
    subl $1, %edx
    xorl %eax, %eax
    leaq 8(,%rdx,8), %rdx
    xorpd %xmm1, %xmm1 /* xmm1 is always 0.0 */
.tl_loopstart:
    movsd (%rsi,%rax), %xmm0  /* xmm0 = pred[i] */
    ucomisd %xmm0, %xmm1
    jae .tl_predzero
    fldl (%rdi,%rax)
    fldl (%rsi,%rax)
    fyl2x
    addsd %xmm0, %xmm2
    faddp %st, %st(1)
.tl_loopinc:
    addq $8, %rax
    cmpq %rdx, %rax
    jne .tl_loopstart
.tl_exit:
    fstpl -8(%rsp)
    movsd -8(%rsp), %xmm0
    mulsd  .minusoneoverlog2e(%rip), %xmm0
    addsd  %xmm2, %xmm0
    ret
.tl_predzero: /* look at data: */
    movsd (%rdi,%rax), %xmm0
    ucomisd %xmm0, %xmm1
    /* if data is also zero, skip this entry and go to the next: */
    je .tl_loopinc
    /* otherwise: return infinity; pop fpu stack first: */
    fstpl -8(%rsp)
    movsd .infinity(%rip), %xmm0
    ret

.section .rodata
.align 8
.infinity:
    .byte 0, 0, 0, 0, 0, 0, 0xf0, 0x7f
.minusoneoverlog2e: /* -1.0 / log_2(e)  */
    .byte 239, 57, 250, 254, 66, 46, 230, 191

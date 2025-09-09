<!-- TestFormula.vue -->
<template>
  <div>
    <h2>行内公式：</h2>
    <div v-html="inlineFormula"></div>

    <h2>块级公式：</h2>
    <div v-html="blockFormula"></div>
  </div>
</template>

<script>
import katex from 'katex'

export default {
  data() {
    return {
      inlineFormula: '',
      blockFormula: ''
    }
  },
  mounted() {
    // 行内公式
    this.inlineFormula = katex.renderToString('c = \\pm\\sqrt{a^2 + b^2}', {
      throwOnError: false,
      displayMode: false // 行内模式
    })

    // 块级公式
    this.blockFormula = katex.renderToString('\\int_{-\\infty}^\\infty e^{-x^2} \\, dx = \\sqrt{\\pi}', {
      throwOnError: false,
      displayMode: true // 块级模式（居中、独立行）
    })
  }
}
</script>

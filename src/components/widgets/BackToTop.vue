<template>
  <button v-show="showBackToTop" class="back-to-top" @click="backToTop" aria-label="回到顶部">
    <i class="fas fa-arrow-up"></i>
  </button>
</template>

<script>
export default {
  name: 'BackToTop',
  data() {
    return {
      showBackToTop: false
    }
  },
  mounted() {
    window.addEventListener('scroll', this.handleScroll)
    this.handleScroll()
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.handleScroll)
  },
  methods: {
    handleScroll() {
      this.showBackToTop = window.scrollY > 180
    },
    backToTop() {
      window.scrollTo({ top: 0, behavior: 'smooth' })
    }
  }
}
</script>

<style scoped>
.back-to-top {
  position: fixed;
  bottom: 30px;
  right: 30px;
  width: 40px;
  height: 40px;
  background: #fff;
  color: #495057;
  border: 1px solid #dee2e6;
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 18px;
  font-weight: bold;
  cursor: pointer;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
  z-index: 1000;
}

.back-to-top:hover {
  background: #f8f9fa;
  border-color: #adb5bd;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.12);
}

.back-to-top:active {
  transform: scale(0.95);
}
</style>

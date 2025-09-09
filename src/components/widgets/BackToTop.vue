<template>
  <button v-show="showBackToTop" class="back-to-top" @click="handleClick" aria-label="回到顶部"
    @touchstart="handleTouchStart" @touchmove="handleTouchMove" @touchend="handleTouchEnd"
    :style="{ top: buttonTop + 'px' }">
    <i class="fas fa-arrow-up"></i>
  </button>
</template>

<script>
export default {
  name: 'BackToTop',
  data() {
    return {
      showBackToTop: false,
      isDragging: false,
      startY: 0,
      initialTop: 0,
      buttonTop: window.innerHeight - 100,
      touchMoved: false
    }
  },
  mounted() {
    window.addEventListener('scroll', this.handleScroll)
    this.handleScroll()
    this.buttonTop = window.innerHeight - 100
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.handleScroll)
  },
  methods: {
    handleScroll() {
      if (!this.isDragging) {
        this.showBackToTop = window.scrollY > 180
      }
    },
    backToTop() {
      window.scrollTo({ top: 0, behavior: 'smooth' })
    },
    handleClick() {
      if (!this.isDragging) {
        this.backToTop()
      }
    },
    handleTouchStart(e) {
      this.isDragging = true
      this.touchMoved = false
      this.startY = e.touches[0].clientY
      this.initialTop = this.buttonTop
    },
    handleTouchMove(e) {
      this.touchMoved = true

      const currentY = e.touches[0].clientY
      const diffY = currentY - this.startY

      let newTop = this.initialTop + diffY

      const buttonHeight = 40
      const minTop = 20
      const maxTop = window.innerHeight - buttonHeight - 20

      newTop = Math.max(minTop, newTop)
      newTop = Math.min(maxTop, newTop)

      this.buttonTop = newTop

      e.preventDefault()
    },
    handleTouchEnd() {
      this.isDragging = false

      if (!this.touchMoved) {
        this.backToTop()
      }
    }
  }
}
</script>

<style scoped>
.back-to-top {
  position: fixed;
  right: 30px;
  width: 40px;
  height: 40px;
  background: #fff;
  color: #495057;
  border: none;
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 18px;
  font-weight: bold;
  cursor: pointer;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
}

.back-to-top:hover {
  background: #f8f9fa;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.12);
}

.back-to-top:active {
  transform: scale(0.95);
}
</style>

<template>
  <button
    v-show="visible"
    class="toc-drawer-btn d-lg-none d-flex align-items-center justify-content-center"
    @click="openDrawer"
    :aria-label="t('openToc')"
    @touchstart.prevent.stop="onTouchStart"
    @touchmove.prevent.stop="onTouchMove"
    @touchend.prevent.stop="onTouchEnd"
    :style="{ top: buttonTop + 'px' }"
  >
    <i class="fas fa-bookmark"></i>
  </button>

  <div v-if="show" class="mobile-offcanvas d-lg-none" @click.self="close">
    <div class="offcanvas-panel offcanvas-right border-start rounded-0 shadow-sm">
      <div class="offcanvas-section">
        <div class="offcanvas-card">
          <OnThisPage
            containerSelector=".markdown-body"
            :levels="[2, 3]"
            :offset="88"
            @navigate="onNavigate"
          />
        </div>
      </div>
    </div>
    <div class="offcanvas-backdrop" @click="close"></div>
  </div>
</template>

<script setup lang="ts">
import { onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import OnThisPage from '@/components/layout/OnThisPage.vue'
import { useFloatingButton } from '@/composables/useFloatingButton'
import { lockScrollPosition, unlockScrollPosition } from '@/utils/scroll'

const { t } = useI18n()

const visible = ref(false)
const show = ref(false)
const lockedScrollY = ref<number | null>(null)

const {
  isDragging,
  buttonTop,
  clampTop,
  dispatchBaseTop,
  onTouchStart,
  onTouchMove,
  onTouchEnd,
  subscribe,
  unsubscribe,
} = useFloatingButton({
  sourceId: 'toc',
  defaultTop: (typeof window !== 'undefined' ? window.innerHeight : 1024) - 160,
  mode: 'stack',
  onRelease: openDrawer,
})

function onResize() {
  const isMobile = window.innerWidth < 992
  visible.value = isMobile
  buttonTop.value = clampTop(buttonTop.value)
}

function openDrawer() {
  if (isDragging.value) return
  lockedScrollY.value = lockScrollPosition()
  show.value = true
}

function close() {
  show.value = false
  unlockScrollPosition(lockedScrollY.value)
}

function onNavigate() {
  close()
}

onMounted(() => {
  window.addEventListener('resize', onResize, { passive: true })
  subscribe()
  onResize()
  dispatchBaseTop()
})

onBeforeUnmount(() => {
  window.removeEventListener('resize', onResize)
  unsubscribe()
  unlockScrollPosition(lockedScrollY.value)
})
</script>

<style scoped>
.toc-drawer-btn {
  position: fixed;
  right: 28px;
  width: 40px;
  height: 40px;
  background-color: var(--surface);
  color: var(--fg-2);
  border: 1px solid var(--line);
  border-radius: 50%;
  font-size: 14px;
  cursor: pointer;
  box-shadow: var(--shadow-soft);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
  transition:
    background-color 0.14s ease,
    box-shadow 0.14s ease,
    transform 0.14s ease,
    color 0.14s ease,
    border-color 0.14s ease;
}

.toc-drawer-btn:hover {
  background-color: var(--primary);
  color: var(--on-primary);
  border-color: transparent;
  box-shadow: var(--shadow-lift);
}

.toc-drawer-btn:active {
  transform: scale(0.93);
}

.mobile-offcanvas {
  position: fixed;
  inset: 0;
  z-index: 1050;
}

.offcanvas-backdrop {
  position: absolute;
  inset: 0;
  background: var(--app-backdrop-bg);
  z-index: 1;
}

.offcanvas-panel {
  position: absolute;
  top: 0;
  bottom: 0;
  width: min(82vw, 300px);
  background: var(--surface);
  border-left: 1px solid var(--line);
  box-shadow: var(--shadow-lift);
  overflow-y: auto;
  -webkit-overflow-scrolling: touch;
  padding: 1rem 0.75rem 1.5rem;
  z-index: 2;
}

.offcanvas-right {
  right: 0;
}

.offcanvas-section + .offcanvas-section {
  margin-top: 0.75rem;
}

.offcanvas-card {
  background-color: var(--surface);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  padding: 0.5rem;
  margin: 0.5rem 0;
}
</style>

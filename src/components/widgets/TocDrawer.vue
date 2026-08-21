<template>
  <FloatingButton
    class="d-lg-none"
    source-id="toc"
    :default-top="defaultTop"
    mode="stack"
    :on-release="openDrawer"
    :ariaLabel="t('openToc')"
    :show="visible"
  >
    <i class="fas fa-bookmark"></i>
  </FloatingButton>

  <div v-if="show" class="mobile-offcanvas d-lg-none" @click.self="close">
    <div class="offcanvas-panel offcanvas-right border-start rounded-0 shadow-sm">
      <div class="offcanvas-section">
        <div class="offcanvas-card">
          <OnThisPage
            containerSelector=".markdown-body"
            :levels="[2, 3]"
            :offset="HEADER_OFFSET"
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
import FloatingButton from '@/components/widgets/FloatingButton.vue'
import { lockScrollPosition, unlockScrollPosition } from '@/utils/scroll'
import { HEADER_OFFSET } from '@/config'

const { t } = useI18n()

const visible = ref(false)
const show = ref(false)
const lockedScrollY = ref<number | null>(null)

const defaultTop = (typeof window !== 'undefined' ? window.innerHeight : 1024) - 160

function onResize() {
  visible.value = window.innerWidth < 992
}

function openDrawer() {
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
  onResize()
})

onBeforeUnmount(() => {
  window.removeEventListener('resize', onResize)
  unlockScrollPosition(lockedScrollY.value)
})
</script>

<style scoped>
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

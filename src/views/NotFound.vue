<template>
  <div class="page-section notfound-view">
    <h1 class="notfound-code">404</h1>
    <p class="notfound-desc">{{ t('pageNotFound') }}</p>
    <router-link :to="homePath" class="notfound-back">
      <i class="fas fa-arrow-left"></i>{{ t('backHome') }}
    </router-link>
  </div>
</template>

<script setup lang="ts">
defineOptions({
  name: 'NotFoundView',
  head() {
    return { title: '404 - zorrooz’s blog' }
  },
})
import { computed } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRoute } from 'vue-router'

const { t } = useI18n()
const route = useRoute()

const homePath = computed(() => {
  const segment = route.params.locale === 'en' ? 'en' : 'zh'
  return `/${segment}/`
})
</script>

<style scoped>
.notfound-view {
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
  padding-top: var(--sp-24);
  padding-bottom: var(--sp-28);
}

.notfound-code {
  font-family: var(--font-serif);
  font-size: clamp(88px, 20vw, 180px);
  font-weight: 700;
  letter-spacing: -0.04em;
  line-height: 1;
  margin: 0;
  color: var(--fg);
}

.notfound-desc {
  font-size: var(--text-lg);
  color: var(--fg-2);
  margin: var(--sp-6) 0 var(--sp-10);
  line-height: 1.65;
}

.notfound-back {
  display: inline-flex;
  align-items: center;
  gap: var(--sp-3);
  padding: 10px 24px;
  border: 1px solid var(--line);
  border-radius: var(--radius-pill);
  font-size: var(--text-base);
  font-weight: 500;
  color: var(--fg-2);
  background: var(--surface);
  text-decoration: none;
  transition:
    color var(--dur-fast) ease,
    background-color var(--dur-fast) ease,
    border-color var(--dur-fast) ease;
}

.notfound-back:hover {
  border-color: color-mix(in srgb, var(--primary) 40%, transparent);
  color: var(--primary);
  background: var(--tint);
}

.notfound-back i {
  font-size: 13px;
}
</style>
